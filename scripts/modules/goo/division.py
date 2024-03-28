import bpy, bmesh
from goo.cell import Cell


class DivisionLogic:
    def make_divide(self, mother: Cell) -> tuple[Cell, Cell]:
        pass


class BisectDivisionLogic(DivisionLogic):
    def make_divide(self, mother):
        com = mother.COM(global_coords=False)
        axis = mother.major_axis().axis()

        daughter = mother.copy()
        self._bisect(mother.obj.data, com, axis, True)
        self._bisect(daughter.obj.data, com, axis, False)

        mother.remesh()
        daughter.remesh()

        daughter.name = mother.name + ".1"
        mother.name = mother.name + ".0"

        return mother, daughter

    def _bisect(self, mesh, com, axis, inner):
        bm = bmesh.new()
        bm.from_mesh(mesh)

        # bisect with plane
        verts = [v for v in bm.verts]
        edges = [e for e in bm.edges]
        faces = [f for f in bm.faces]
        geom = verts + edges + faces

        result = bmesh.ops.bisect_plane(
            bm,
            geom=geom,
            plane_co=com,
            plane_no=axis,
            clear_outer=inner,
            clear_inner=not inner,
        )

        # fill in bisected face
        edges = [e for e in result["geom_cut"] if isinstance(e, bmesh.types.BMEdge)]
        bmesh.ops.edgeloop_fill(bm, edges=edges)

        bm.to_mesh(mesh)
        bm.free()
        mesh.update()


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
        # TODO: this is expensive, and requires disabling physics/evaluating depsgraph for each call. Maybe look at different contexts, or creating new objects?
        bpy.ops.object.modifier_apply(modifier=bool_mod.name)

        # separate two daughter cells
        # TODO: ops are expensive, look to reduce this to low-level.
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
