import bpy, bmesh
from goo.cell import Cell


class DividerHandler:
    def divide(self, mother: Cell) -> tuple[Cell, Cell]:
        pass


class BooleanDivider(DividerHandler):
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
        # TODO: ops are expensive, look to reduce this to low-level
        bpy.ops.object.mode_set(mode="EDIT")
        bpy.ops.mesh.separate(type="LOOSE")
        bpy.ops.object.mode_set(mode="OBJECT")

        daughter = Cell(bpy.context.selected_objects[0])
        daughter.name = mother.name + ".1"
        mother.name = mother.name + ".0"

        # remesh daughter cells
        self._division_remesh(mother.obj)
        self._division_remesh(daughter.obj)

        # clean up
        mesh = plane.data
        bpy.data.meshes.remove(mesh, do_unlink=True)

        return mother, daughter

    def _division_remesh(self, obj):
        bpy.context.view_layer.objects.active = obj
        remesh_modifier = obj.modifiers.new(name="Remesh", type="REMESH")
        remesh_modifier.mode = "VOXEL"
        remesh_modifier.voxel_size = 0.25  # microns
        remesh_modifier.adaptivity = 0
        remesh_modifier.use_remove_disconnected = True
        # remesh_modifier.use_smooth_shade = True
        bpy.ops.object.modifier_apply(modifier="Remesh")
