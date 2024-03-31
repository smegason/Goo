import bpy, bmesh
from goo.cell import Cell
import numpy as np


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
        com = mother.COM(global_coords=False)
        axis = mother.major_axis().axis()

        daughter = mother.copy()

        m_mb = self._bisect(mother.obj_eval, com, axis, True, self.margin)
        d_mb = self._bisect(mother.obj_eval.copy(), com, axis, False, self.margin)
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
        self.to_flush.clear()


class BooleanDivisionLogic(DivisionLogic):
    def __init__(self):
        pass

    def make_divide(self, mother):
        plane = mother.create_division_plane()
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
