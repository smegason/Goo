bl_info = {
    "name": "Goo - Blender add-on",
    "author": "Antoine Ruzette",
    "version": (0, 1),
    "blender": (3, 2, 0),
    "location": "View3D > N",
    "description": "Creates user-specified cell-like objects",
    "warning": "",
    "doc_url": "https://smegason.github.io/Goo/",
    "category": "Add Random",
}

import bpy
from goo import goo_the_blender_way

def create_cell(context):

    cA_collection = bpy.data.collections.new("A_Cells")
    bpy.context.scene.collection.children.link(cA_collection)
    goo_the_blender_way.make_cell(name='CellA1', location=(0,0,0), radius=1, size=1, scale=(1,1,1))
    obj = bpy.context.active_object
    bpy.ops.collection.objects_remove_all()
    # Add the active cell to our specific collection 
    bpy.data.collections['A_Cells'].objects.link(obj)


class AddCellOperator(bpy.types.Operator):
    """Tooltip"""
    bl_idname = "create.cell"
    bl_label = "Create cell"

    def execute(self, context):
        create_cell(context)
        return {'FINISHED'}


class GooPanel(bpy.types.Panel):
    """Creates a Panel in the scene context of the properties editor"""
    bl_label = "Goo"
    bl_idname = "SCENE_PT_layout"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "scene"

    def draw(self, context):
        layout = self.layout
        scene = context.scene

        '''        
        # Create a simple row.
        layout.label(text=" Simple Row:")

        row = layout.row()
        row.prop(scene, "frame_start")
        row.prop(scene, "frame_end")

        # Create an row where the buttons are aligned to each other.
        layout.label(text=" Aligned Row:")

        row = layout.row(align=True)
        row.prop(scene, "frame_start")
        row.prop(scene, "frame_end")

        # Create two columns, by using a split layout.
        split = layout.split()

        # First column
        col = split.column()
        col.label(text="Column One:")
        col.prop(scene, "frame_end")
        col.prop(scene, "frame_start")

        # Second column, aligned
        col = split.column(align=True)
        col.label(text="Column Two:")
        col.prop(scene, "frame_start")
        col.prop(scene, "frame_end")
        '''

        # Big render button
        layout.label(text="Create cell:")
        row = layout.row()
        row.scale_y = 3.0
        row.operator("create.cell")


_classes = {
    AddCellOperator,
    GooPanel
}


def register():
    for cls in _classes:
        bpy.utils.register_class(cls)

def unregister():
    for cls in _classes:
        bpy.utils.unregister_class(cls)

if __name__ == "__main__":
    register()
