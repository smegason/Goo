from typing import Callable

import bpy, bmesh
from mathutils import Vector
import numpy as np

from goo.cell import Cell
from goo.utils import *
from goo.handler import Handler


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
        com = mother.COM(local_coords=True)
        axis = mother.major_axis().axis(local_coords=True)
        base_name = mother.name

        m_mb = self._bisect(mother.obj_eval, com, axis, True, self.margin)
        d_mb = self._bisect(mother.obj_eval, com, axis, False, self.margin)

        daughter = mother.copy()

        mother.name = base_name + ".0"
        daughter.name = base_name + ".1"

        self.to_flush.append((m_mb, mother))
        self.to_flush.append((d_mb, daughter))

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

        plane = self._create_division_plane(
            mother.name, mother.major_axis(), mother.COM()
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

    def _create_division_plane(self, name, major_axis, com, collection=None):
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
        bpy.context.scene.collection.objects.link(plane)

        # Add thickness to plane
        solid_mod = plane.modifiers.new(name="Solidify", type="SOLIDIFY")
        solid_mod.offset = 0
        solid_mod.thickness = 0.025

        plane.hide_set(True)
        return plane


class DivisionHandler(Handler):
    def __init__(self, divider_logic):
        self.divider_logic = divider_logic()

    def setup(self, get_cells: Callable[[], list[Cell]], dt):
        super(DivisionHandler, self).setup(get_cells, dt)
        for cell in self.get_cells():
            cell["divided"] = False
        self._cells_to_update = []

    def can_divide(self, cell: Cell) -> bool:
        raise NotImplementedError("Subclasses must implement can_divide() method.")

    def update_on_divide(self, cell: Cell):
        pass

    def run(self, scene, depsgraph):
        for cell in self._cells_to_update:
            cell.enable_physics()
            cell.cloth_mod.point_cache.frame_start = scene.frame_current
            cell["divided"] = False
        self._cells_to_update.clear()

        for cell in self.get_cells():
            if self.can_divide(cell):
                mother, daughter = cell.divide(self.divider_logic)
                self.update_on_divide(mother)
                self.update_on_divide(daughter)

                if mother.physics_enabled:
                    self._cells_to_update.extend([mother, daughter])

        for cell in self._cells_to_update:
            cell.disable_physics()
            cell["divided"] = True
        self.divider_logic.flush()


class TimeDivisionHandler(DivisionHandler):
    def __init__(self, divider_logic, mu=10, var=0):
        super(TimeDivisionHandler, self).__init__(divider_logic)
        self.mu = mu
        self.var = var

    def setup(self, get_cells, dt):
        super(TimeDivisionHandler, self).setup(get_cells, dt)
        for cell in self.get_cells():
            cell["last_division_time"] = 0

    # TODO: implement variance
    def can_divide(self, cell: Cell):
        time = bpy.context.scene.frame_current * self.dt
        if "last_division_time" not in cell:
            cell["last_division_time"] = time
            return False
        return time - cell["last_division_time"] >= self.mu

    def update_on_divide(self, cell: Cell):
        time = bpy.context.scene.frame_current * self.dt
        cell["last_division_time"] = time


class SizeDivisionHandler(DivisionHandler):
    def __init__(self, divider_logic, threshold=30):
        super(SizeDivisionHandler, self).__init__(divider_logic)
        self.threshold = threshold

    def can_divide(self, cell: Cell):
        return cell.volume() >= self.threshold
