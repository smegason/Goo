import bpy, bmesh
from mathutils import Vector
import numpy as np

from goo.cell import Cell
from goo.utils import *


class DivisionLogic:
    def make_divide(self, mother: Cell) -> tuple[Cell, Cell]:
        pass

    def flush(self):
        pass


class BisectDivisionLogic(DivisionLogic):
    def __init__(self, margin=0.025):
        self.margin = margin
        self.to_flush = []

    def make_divide(self, mother):
        com = mother.get_COM(local_coords=True)
        axis = mother.get_major_axis().axis(local_coords=True)
        obj_eval = mother.obj_eval.copy()

        daughter = mother.copy()

        m_mb = self._bisect(obj_eval, com, axis, True, self.margin)
        d_mb = self._bisect(obj_eval, com, axis, False, self.margin)

        self.to_flush.append((m_mb, mother))
        self.to_flush.append((d_mb, daughter))

        daughter.name = mother.name + ".1"
        mother.name = mother.name + ".0"

        return mother, daughter

    def _bisect(self, obj_eval, com, axis, inner, margin):
        bm = bmesh.new()
        bm.from_mesh(obj_eval.to_mesh())

        # bisect with plane
        verts = [v for v in bm.verts]
        edges = [e for e in bm.edges]
        faces = [f for f in bm.faces]
        geom = verts + edges + faces

        result = bmesh.ops.bisect_plane(
            bm,
            geom=geom,
            plane_co=com + axis * margin / 2 if inner else com - axis * margin / 2,
            plane_no=axis,
            clear_inner=inner,
            clear_outer=not inner,
        )

        # fill in bisected face
        edges = [e for e in result["geom_cut"] if isinstance(e, bmesh.types.BMEdge)]
        bmesh.ops.edgeloop_fill(bm, edges=edges)

        return bm

    def flush(self):
        for bm, cell in self.to_flush:
            bm.to_mesh(cell.obj.data)
            bm.free()
            cell.remesh()
            cell.recenter()
        self.to_flush.clear()


class BooleanDivisionLogic(DivisionLogic):
    def __init__(self):
        pass

    def make_divide(self, mother: Cell):

        plane = self.create_division_plane(
            mother.name, mother.get_major_axis(), mother.get_COM()
        )
        obj = mother.obj

        # cut mother cell by division plane
        bpy.context.view_layer.objects.active = obj
        bool_mod = obj.modifiers.new(name="Boolean", type="BOOLEAN")
        bool_mod.operand_type = "OBJECT"
        bool_mod.object = plane
        bool_mod.operation = "DIFFERENCE"
        bool_mod.solver = "EXACT"
        bpy.ops.object.modifier_apply(modifier=bool_mod.name)

        # separate two daughter cells
        bpy.ops.object.mode_set(mode="EDIT")
        bpy.ops.mesh.separate(type="LOOSE")
        bpy.ops.object.mode_set(mode="OBJECT")

        daughter = Cell(bpy.context.selected_objects[0])
        daughter.obj.select_set(False)

        daughter.name = mother.name + ".1"
        mother.name = mother.name + ".0"

        # remesh daughter cells
        mother.remesh()
        daughter.remesh()

        # clean up
        bpy.data.meshes.remove(plane.data, do_unlink=True)

        return mother, daughter

    def flush(self):
        pass

    def _create_division_plane(name, major_axis, com, collection=None):
        """
        Creates a plane orthogonal to the long axis vector
        and passing through the cell's center of mass.
        """
        # Define new plane
        plane = create_mesh(
            f"{name}_division_plane",
            loc=com,
            mesh="plane",
            size=major_axis.length() + 1,
            rotation=major_axis.axis().to_track_quat("Z", "Y"),
        )

        # Add thickness to plane
        solid_mod = plane.modifiers.new(name="Solidify", type="SOLIDIFY")
        solid_mod.offset = 0
        solid_mod.thickness = 0.025

        plane.hide_set(True)
        return plane


class TimeDivisionHandler:
    # TODO: implement variance
    def __init__(self, divider_handler, mu=10, var=0):
        self.mu = mu
        self.var = var
        self.divider_handler = divider_handler()

    def setup(self, get_cells, dt):
        self.get_cells = get_cells
        self.dt = dt

    def run(self, scene, depsgraph):
        time = scene.frame_current * self.dt
        cells = self.get_cells()
        for cell in self.get_cells():
            if time - cell.last_division_time >= self.mu:
                mother, daughter = cell.divide(self.divider_handler)
                mother.last_division_time = time
                daughter.last_division_time = time
        self.divider_handler.flush()


class TimeDivisionPhysicsHandler(TimeDivisionHandler):
    def setup(self, get_cells, dt):
        super(TimeDivisionPhysicsHandler, self).setup(get_cells, dt)
        self._cells_to_update = []

    def run(self, scene, depsgraph):
        for cell in self._cells_to_update:
            cell.enable_physics(collision=False)
            cell.cloth_mod.point_cache.frame_start = scene.frame_current
        self._cells_to_update.clear()

        time = scene.frame_current * self.dt
        for cell in self.get_cells():
            if time - cell.last_division_time >= self.mu:
                mother, daughter = cell.divide(self.divider_handler)
                self._cells_to_update.extend([mother, daughter])
        for cell in self._cells_to_update:
            cell.disable_physics(collision=False)
            cell.last_division_time = time
        self.divider_handler.flush()
