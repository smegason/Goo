import bpy
from goo import goo 
from importlib import reload
reload(goo)

# Delete all existing collections 
for collection in bpy.data.collections:  # loop through the existing collection
    # loop through objects in collection
    for objs in collection.objects:
        # delete existing objects in collection 
        bpy.data.objects.remove(objs)
    # Delete collection
    bpy.data.collections.remove(collection)

# Change the Viewport Shading to Rendered
for area in bpy.data.screens[3].areas:
    if area.type == 'VIEW_3D':
        for space in area.spaces: 
            if space.type == 'VIEW_3D':
                space.shading.type = 'MATERIAL'  # 'RENDERED', 'MATERIAL'


# ====================================================
cA_collection = bpy.data.collections.new("A_Cells")
# link the collection to the scene for visualization 
bpy.context.scene.collection.children.link(cA_collection)

# Add material 
goo.add_material_cell("red", 1,0,0)

# Define cell A1
cA1 = goo.Cell("cell_A1", loc=(0, 0, 0), material="red", flavor="ico_sphere")
# Make a Blender mesh object for cell
goo.make_cell(cA1)

# Calculate volume of an object
obj = bpy.data.objects[cA1.data['name']]
vol = goo.calculate_volume(obj)
#print(vol)

# Calculates the major axis 
axis =  goo.get_major_axis(obj)
#print(axis)

# Moves cell along the given axis 
axis = (20, 20, 0)
goo.translate_cell(obj, axis)

# The created cell is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene 
# if you do not delete them the object outside colleciton will link to object  
# inside the collection
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection
bpy.data.collections['A_Cells'].objects.link(obj)
