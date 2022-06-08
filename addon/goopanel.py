bl_info = {
    "name": "Custom addon",
    "author": "Sean Megason",
    "version": (1, 0),
    "blender": (2, 80, 0),
    "location": "View3D > N",
    "description": "Adds random spheres",
    "warning": "",
    "doc_url": "",
    "category": "Add Random",
}


import bpy
from goo import goo
from bpy.types import (Panel, Operator)

class AddCellOperator(bpy.types.Operator):
    """Constructs a cell at the cursor"""
    bl_idname = "addcell.1"
    bl_label = "Add Cell"

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

from bpy.utils import register_class, unregister_class

_classes = {
    AddCellOperator,
    CustomPanel
}

def register():
    for cls in _classes:
        register_class(cls)

def unregister():
    for cls in _classes:
        unregister_class(cls)

if __name__ == "__main__":
    register()
