import bpy
from goo import goo
from mathutils import *
D = bpy.data
C = bpy.context

# Delete all existing collections 
for collection in bpy.data.collections:  # loop through the existing collection
    # loop through objects in collection
    for objs in collection.objects:
        # delete existing objects in collection 
        bpy.data.objects.remove(objs)
    # Delete collection
    bpy.data.collections.remove(collection)

cA_collection = bpy.data.collections.new("A_Cells")
bpy.context.scene.collection.children.link(cA_collection)

cell = goo.Cell("cell_A1", loc=(2, 2, 0))
goo.make_cell(cell)

obj = bpy.context.active_object
bpy.ops.collection.objects_remove_all()
bpy.data.collections['A_Cells'].objects.link(obj)