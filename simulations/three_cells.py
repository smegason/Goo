# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy

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
for area in bpy.context.screen.areas: 
    if area.type == 'VIEW_3D':
        for space in area.spaces: 
            if space.type == 'VIEW_3D':
                space.shading.type = 'MATERIAL'
                
#====================== Colors =========================
# Green
matg = bpy.data.materials.new("Green")
matg.diffuse_color = (0,0.1,0,0.8)
# Red
matr = bpy.data.materials.new("red")
matr.diffuse_color = (0.1,0,0,0.8)
# Blue
matb = bpy.data.materials.new("blue")
matb.diffuse_color = (0.0,0,1,0.8)

#================== Cell A Collection ==================
# Create a collection for cell A
cA_collection = bpy.data.collections.new("A_Cells")
# link the collection to the scene for visualization 
bpy.context.scene.collection.children.link(cA_collection)
# Define cell A1
cA1 = goo.Cell("cell_A1", loc = (2,2,0))
# Make a Blender mesh object for cell
goo.make_cell(cA1)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene 
# if you do not delete them the object outside colleciton will link to object inside the collection
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)


#================== Cell B Collection ==================
# Create a collection for cell B
cB_collection = bpy.data.collections.new("B_Cells")
# link the collection to the scene
bpy.context.scene.collection.children.link(cB_collection)
# Define cell A1
cB1 = goo.Cell("cell_B1", loc = (-2,2,0))
# Make a Blender mesh object for cell
goo.make_cell(cB1)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matr
# Remove object from all collections not used in a scene 
# if you do not delete them the object outside colleciton will link to object inside the collection
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['B_Cells'].objects.link(obj)

#================== Cell C Collection ==================
# Create a collection for cell B
cC_collection = bpy.data.collections.new("C_Cells")
# link the collection to the scene
bpy.context.scene.collection.children.link(cC_collection)
# Define cell A1
cC1 = goo.Cell("cell_C1", loc = (0,2,2))
# Make a Blender mesh object for cell
goo.make_cell(cC1)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matb
# Remove object from all collections not used in a scene 
# if you do not delete them the object outside colleciton will link to object inside the collection
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['C_Cells'].objects.link(obj)

#================== Force A Collection ==================
# Create a collection for force A
fA_collection = bpy.data.collections.new("A_Forces")
# link the collection to the scene 
bpy.context.scene.collection.children.link(fA_collection)
# Define and link force A1
fA1 = goo.Force("force_A1", "cell_A1", -800)
# Make force
goo.make_force(fA1)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

#================== Force B Collection ==================
# Create a collection for force B
fB_collection = bpy.data.collections.new("B_Forces")
# link the collection to the scene 
bpy.context.scene.collection.children.link(fB_collection)
# Define and link force B1
fB1 = goo.Force("force_B1", "cell_B1", -800)
# Make force
goo.make_force(fB1)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['B_Forces'].objects.link(obj)

#================== Force C Collection ==================
# Create a collection for force B
fC_collection = bpy.data.collections.new("C_Forces")
# link the collection to the scene 
bpy.context.scene.collection.children.link(fC_collection)
# Define and link force B1
fC1 = goo.Force("force_C1", "cell_C1", -800)
# Make force
goo.make_force(fC1)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['C_Forces'].objects.link(obj)


# Set the Effector Collection (Select a cell then go to
# Physics Properties > Cloth > Field Weights > Effector Collection > selct the force collection)
for collection in bpy.data.collections:
    
    # Exclude the objects in the force collections
    if 'Cells' in collection.name_full:
        # Collection name
        coll_name = collection.name_full
        
        # Cell A attracts to Cell B
        if coll_name == 'A_Cells':
                      
            # Force B
            force_name = 'B_Forces'
            # Loop through the objects existed in the collection 
            for obj in collection.objects:
                obj.modifiers['Cloth'].settings.effector_weights.collection = bpy.data.collections[force_name]
            
        # Cell B attracts to Cell C
        elif coll_name == 'B_Cells':
            
            # Force C
            force_name = 'C_Forces'
            # Loop through the objects existed in the collection 
            for obj in collection.objects:
                obj.modifiers['Cloth'].settings.effector_weights.collection = bpy.data.collections[force_name]
        
        # Cell C attracts to Cell A    
        else:
            
            # Foce A
            force_name = 'A_Forces'
            # Loop through the objects existed in the collection 
            for obj in collection.objects:
                obj.modifiers['Cloth'].settings.effector_weights.collection = bpy.data.collections[force_name]
                

        
handlers = goo.handler_class()
handlers.forces = [fA1, fB1, fC1]
bpy.app.handlers.frame_change_post.clear()
bpy.app.handlers.frame_change_post.append(handlers.adhesion_handler)
