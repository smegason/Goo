bl_info = {
    "name": "Goo - Blender add-on",
    "author": "Antoine Ruzette, Sean Megason",
    "version": (0, 1),
    "blender": (3, 2, 0),
    "location": "View3D > N",
    "description": "Creates user-specified cell-like objects",
    "warning": "",
    "doc_url": "",
    "category": "Add Random",
}


import bpy
from goo import goo

class AddCellOperator(bpy.types.Operator):
    """Constructs a cell at the cursor"""
    bl_idname = "addcell.1"
    bl_label = "Add Cell"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        goo.setup_world()
        center = bpy.context.scene.cursor.location
        cell = goo.Cell(name_string="Cell1_", loc=center, material="CellGreen")
        goo.make_cell(cell)

        return {'FINISHED'}

class CustomPanel(bpy.types.Panel):
    """Creates a Panel in the Object properties window"""
    bl_label = "Goo"
    bl_idname = "OBJECT_PT_random"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "Goo"

    def draw(self, context):
        layout = self.layout
        obj = context.object
        row = layout.row()
        row.operator(AddCellOperator.bl_idname, text=AddCellOperator.bl_label, icon='SPHERE')

_classes = {
    AddCellOperator,
    CustomPanel
}

def register():
    for cls in _classes:
        bpy.utils.register_class(cls)

def unregister():
    for cls in _classes:
        bpy.utils.unregister_class(cls)

if __name__ == "__main__":
    register()
