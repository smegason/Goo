import bpy
from importlib import reload
from goo import goo_the_blender_way
from goo import goo

reload(goo_the_blender_way)
reload(goo)
goo.setup_world()
    
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
                space.shading.type = 'RENDERED'


#================== Cell Green Collection ==================
# Create a collection for cell Green
green_collection = goo_the_blender_way.create_collection("Green_cells")
# Define cell G1
g1 = goo_the_blender_way.make_cell_blender(name = 'cell_G1', radius = 1, location = (1,1,1), material = 'Green')
# Add the active cell to our specific collection 
goo_the_blender_way.link_obj2collection(g1, green_collection)

# Define cell G2
g2 = goo_the_blender_way.make_cell_blender(name = 'cell_G2', radius = 1, location = (-1,1,1), material = 'Green')
# Add the active cell to our specific collection 
goo_the_blender_way.link_obj2collection(g2, green_collection)

# radius and vertices are custom properties
print(f'----------------------------- \n \
        {g1.name}; \n \
        {g1.active_material}; \n \
        {g1.location}; \n \
        {g1.users_collection}; \n \
        {g1["radius"]}; \n \
        {g1["vertices_local"]}; \n \
        {g1["vertices_global"]}; \n \
        {g1["vertices_coord_test"]}; '
        )

'''
# Define cell G2
goo_the_blender_way.make_cell_blender(name = 'cell_G1', radius = 1, location = (1,1,1), material = 'Green')
# The created cell is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene 
# if you do not delete them the object outside colleciton will link to object inside the collection
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['Green_cells'].objects.link(obj)
'''

#================== Cell Red Collection ==================
# Create a collection for cell R1
red_collection = goo_the_blender_way.create_collection("Red_cells")
# Define cell R1
r1 = goo_the_blender_way.make_cell_blender(name = 'cell_R1', radius = 1, location = (3,3,1), material = 'red')
# Add the active cell to our specific collection 
goo_the_blender_way.link_obj2collection(r1, red_collection)

# Define cell R2
r2 = goo_the_blender_way.make_cell_blender(name = 'cell_R2', radius = 1, location = (1,3,1), material = 'red')
# Add the active cell to our specific collection 
goo_the_blender_way.link_obj2collection(r2, red_collection)


'''
bpy.app.handlers.frame_change_post.clear()
bpy.app.handlers.frame_change_post.append(goo_the_blender_way.data_handler)
'''