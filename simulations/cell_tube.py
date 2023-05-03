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
for area in bpy.data.screens[3].areas: 
    if area.type == 'VIEW_3D':
        for space in area.spaces: 
            if space.type == 'VIEW_3D':
                space.shading.type = 'RENDERED'

#====================== Colors =========================
# Green
matg = bpy.data.materials.new("Green")
matg.diffuse_color = (0,0.1,0,0.8)
# Red
matr = bpy.data.materials.new("red")
matr.diffuse_color = (0.1,0,0,0.8)

force_strength = -1000
force_falloff = 1


#================== Cell A Collection ==================
# Create a collection for cell A
cA_collection = bpy.data.collections.new("A_Cells")
# link the collection to the scene for visualization 
bpy.context.scene.collection.children.link(cA_collection)
# Define cell A11
cA11 = goo.Cell("cell_A11", loc = (0,0,3))
# Make a Blender mesh object for cell
goo.make_cell(cA11)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene 
# if you do not delete them the object outside colleciton will link to object inside the collection
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A12
cA12 = goo.Cell("cell_A12", loc = (2.2,0,2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA12)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A3
cA13 = goo.Cell("cell_A13", loc = (3,0,0))
# Make a Blender mesh object for cell
goo.make_cell(cA13)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A14
cA14 = goo.Cell("cell_A14", loc = (2.2,0,-2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA14)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A15
cA15 = goo.Cell("cell_A15", loc = (0,0,-3))
# Make a Blender mesh object for cell
goo.make_cell(cA15)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A16
cA16 = goo.Cell("cell_A16", loc = (-2.2,0,-2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA16)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A17
cA17 = goo.Cell("cell_A17", loc = (-3,0,0))
# Make a Blender mesh object for cell
goo.make_cell(cA17)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A18
cA18 = goo.Cell("cell_A18", loc = (-2.2,0,2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA18)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)


# Define cell A21
cA21 = goo.Cell("cell_A21", loc = (0,1.5,3))
# Make a Blender mesh object for cell
goo.make_cell(cA21)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene 
# if you do not delete them the object outside colleciton will link to object inside the collection
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A22
cA22 = goo.Cell("cell_A22", loc = (2.2,1.5,2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA22)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A23
cA23 = goo.Cell("cell_A23", loc = (3,1.5,0))
# Make a Blender mesh object for cell
goo.make_cell(cA23)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene 
# if you do not delete them the object outside colleciton will link to object inside the collection
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A24
cA24 = goo.Cell("cell_A24", loc = (2.2,1.5,-2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA24)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A25
cA25 = goo.Cell("cell_A25", loc = (0,1.5,-3))
# Make a Blender mesh object for cell
goo.make_cell(cA25)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A14
cA26 = goo.Cell("cell_A26", loc = (-2.2,1.5,-2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA26)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A27
cA27 = goo.Cell("cell_A27", loc = (-3,1.5,0))
# Make a Blender mesh object for cell
goo.make_cell(cA27)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A28
cA28 = goo.Cell("cell_A28", loc = (-2.2,1.5,2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA28)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)



# Define cell A31
cA31 = goo.Cell("cell_A31", loc = (0,3,3))
# Make a Blender mesh object for cell
goo.make_cell(cA31)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene 
# if you do not delete them the object outside colleciton will link to object inside the collection
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A32
cA32 = goo.Cell("cell_A32", loc = (2.2,3,2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA32)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A33
cA33 = goo.Cell("cell_A33", loc = (3,3,0))
# Make a Blender mesh object for cell
goo.make_cell(cA33)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene 
# if you do not delete them the object outside colleciton will link to object inside the collection
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A34
cA34 = goo.Cell("cell_A34", loc = (2.2,3,-2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA34)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A35
cA35 = goo.Cell("cell_A35", loc = (0,3,-3))
# Make a Blender mesh object for cell
goo.make_cell(cA35)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A36
cA36 = goo.Cell("cell_A36", loc = (-2.2,3,-2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA36)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A37
cA37 = goo.Cell("cell_A37", loc = (-3,3,0))
# Make a Blender mesh object for cell
goo.make_cell(cA37)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A38
cA38 = goo.Cell("cell_A38", loc = (-2.2,3,2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA38)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)




# Define cell A41
cA41 = goo.Cell("cell_A41", loc = (0,4.5,3))
# Make a Blender mesh object for cell
goo.make_cell(cA41)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene 
# if you do not delete them the object outside colleciton will link to object inside the collection
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A42
cA42 = goo.Cell("cell_A42", loc = (2.2,4.5,2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA42)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A43
cA43 = goo.Cell("cell_A43", loc = (3,4.5,0))
# Make a Blender mesh object for cell
goo.make_cell(cA43)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene 
# if you do not delete them the object outside colleciton will link to object inside the collection
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A44
cA44 = goo.Cell("cell_A44", loc = (2.2,4.5,-2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA44)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A45
cA45 = goo.Cell("cell_A45", loc = (0,4.5,-3))
# Make a Blender mesh object for cell
goo.make_cell(cA45)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A46
cA46 = goo.Cell("cell_A46", loc = (-2.2,4.5,-2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA46)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A47
cA47 = goo.Cell("cell_A47", loc = (-3,4.5,0))
# Make a Blender mesh object for cell
goo.make_cell(cA47)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A48
cA48 = goo.Cell("cell_A48", loc = (-2.2,4.5,2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA48)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)




# Define cell A51
cA51 = goo.Cell("cell_A51", loc = (0,6,3))
# Make a Blender mesh object for cell
goo.make_cell(cA51)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene 
# if you do not delete them the object outside colleciton will link to object inside the collection
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A52
cA52 = goo.Cell("cell_A52", loc = (2.2,6,2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA52)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A53
cA53 = goo.Cell("cell_A53", loc = (3,6,0))
# Make a Blender mesh object for cell
goo.make_cell(cA53)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A54
cA54 = goo.Cell("cell_A54", loc = (2.2,6,-2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA54)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A55
cA55 = goo.Cell("cell_A55", loc = (0,6,-3))
# Make a Blender mesh object for cell
goo.make_cell(cA55)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A56
cA56 = goo.Cell("cell_A56", loc = (-2.2,6,-2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA56)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A57
cA57 = goo.Cell("cell_A57", loc = (-3,6,0))
# Make a Blender mesh object for cell
goo.make_cell(cA57)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A58
cA58 = goo.Cell("cell_A58", loc = (-2.2,6,2.2))
# Make a Blender mesh object for cell
goo.make_cell(cA58)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)


#================== Force A Collection ==================
# Create a collection for force A
fA_collection = bpy.data.collections.new("A_Forces")
# link the collection to the scene 
bpy.context.scene.collection.children.link(fA_collection)
# Define and link force A11
fA11 = goo.Force("force_A11", "cell_A11", force_strength, force_falloff)
# Make force
goo.make_force(fA11)
# The created force is the active object
obj = bpy.context.active_object

# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A12
fA12 = goo.Force("force_A12", "cell_A12", force_strength, force_falloff)
# Make force
goo.make_force(fA12)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A13
fA13 = goo.Force("force_A13", "cell_A13", force_strength, force_falloff)
# Make force
goo.make_force(fA13)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A14
fA14 = goo.Force("force_A14", "cell_A14", force_strength, force_falloff)
# Make force
goo.make_force(fA14)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A15
fA15 = goo.Force("force_A15", "cell_A15", force_strength, force_falloff)
# Make force
goo.make_force(fA15)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A16
fA16 = goo.Force("force_A16", "cell_A16", force_strength, force_falloff)
# Make force
goo.make_force(fA16)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A17
fA17 = goo.Force("force_A17", "cell_A17", force_strength, force_falloff)
# Make force
goo.make_force(fA17)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A18
fA18 = goo.Force("force_A18", "cell_A18", force_strength, force_falloff)
# Make force
goo.make_force(fA18)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A21
fA21 = goo.Force("force_A21", "cell_A21", force_strength, force_falloff)
# Make force
goo.make_force(fA21)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A22
fA22 = goo.Force("force_A22", "cell_A22", force_strength, force_falloff)
# Make force
goo.make_force(fA22)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A23
fA23 = goo.Force("force_A23", "cell_A23", force_strength, force_falloff)
# Make force
goo.make_force(fA23)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A24
fA24 = goo.Force("force_A24", "cell_A24", force_strength, force_falloff)
# Make force
goo.make_force(fA24)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A25
fA25 = goo.Force("force_A25", "cell_A25", force_strength, force_falloff)
# Make force
goo.make_force(fA25)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A26
fA26 = goo.Force("force_A26", "cell_A26", force_strength, force_falloff)
# Make force
goo.make_force(fA26)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A27
fA27 = goo.Force("force_A27", "cell_A27", force_strength, force_falloff)
# Make force
goo.make_force(fA27)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A28
fA28 = goo.Force("force_A28", "cell_A28", force_strength, force_falloff)
# Make force
goo.make_force(fA28)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)




# Define and link force A31
fA31 = goo.Force("force_A31", "cell_A31", force_strength, force_falloff)
# Make force
goo.make_force(fA31)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A32
fA32 = goo.Force("force_A32", "cell_A32", force_strength, force_falloff)
# Make force
goo.make_force(fA32)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A33
fA33 = goo.Force("force_A33", "cell_A33", force_strength, force_falloff)
# Make force
goo.make_force(fA33)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A34
fA34 = goo.Force("force_A34", "cell_A34", force_strength, force_falloff)
# Make force
goo.make_force(fA34)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A35
fA35 = goo.Force("force_A35", "cell_A35", force_strength, force_falloff)
# Make force
goo.make_force(fA35)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A36
fA36 = goo.Force("force_A36", "cell_A36", force_strength, force_falloff)
# Make force
goo.make_force(fA36)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A37
fA37 = goo.Force("force_A37", "cell_A37", force_strength, force_falloff)
# Make force
goo.make_force(fA37)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A38
fA38 = goo.Force("force_A38", "cell_A38", force_strength, force_falloff)
# Make force
goo.make_force(fA38)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)



# Define and link force A41
fA41 = goo.Force("force_A41", "cell_A41", force_strength, force_falloff)
# Make force
goo.make_force(fA41)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A32
fA42 = goo.Force("force_A42", "cell_A42", force_strength, force_falloff)
# Make force
goo.make_force(fA42)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A43
fA43 = goo.Force("force_A43", "cell_A43", force_strength, force_falloff)
# Make force
goo.make_force(fA43)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A44
fA44 = goo.Force("force_A44", "cell_A44", force_strength, force_falloff)
# Make force
goo.make_force(fA44)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A45
fA45 = goo.Force("force_A45", "cell_A45", force_strength, force_falloff)
# Make force
goo.make_force(fA45)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A36
fA46 = goo.Force("force_A46", "cell_A46", force_strength, force_falloff)
# Make force
goo.make_force(fA46)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A37
fA47 = goo.Force("force_A47", "cell_A47", force_strength, force_falloff)
# Make force
goo.make_force(fA47)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A48
fA48 = goo.Force("force_A48", "cell_A48", force_strength, force_falloff)
# Make force
goo.make_force(fA48)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)





# Define and link force A51
fA51 = goo.Force("force_A51", "cell_A51", force_strength, force_falloff)
# Make force
goo.make_force(fA51)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A52
fA52 = goo.Force("force_A52", "cell_A52", force_strength, force_falloff)
# Make force
goo.make_force(fA52)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A53
fA53 = goo.Force("force_A53", "cell_A53", force_strength, force_falloff)
# Make force
goo.make_force(fA53)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A54
fA54 = goo.Force("force_A54", "cell_A54", force_strength, force_falloff)
# Make force
goo.make_force(fA54)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A55
fA55 = goo.Force("force_A55", "cell_A55", force_strength, force_falloff)
# Make force
goo.make_force(fA55)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A56
fA56 = goo.Force("force_A56", "cell_A56", force_strength, force_falloff)
# Make force
goo.make_force(fA56)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A57
fA57 = goo.Force("force_A57", "cell_A57", force_strength, force_falloff)
# Make force
goo.make_force(fA57)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A58
fA58 = goo.Force("force_A58", "cell_A58", force_strength, force_falloff)
# Make force
goo.make_force(fA58)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Set the Effector Collection (Select a cell then go to
# Physics Properties > Cloth > Field Weights > Effector Collection > selct the force collection)
for collection in bpy.data.collections:
    
    # Exclude the objects in the force collections
    if 'Cells' in collection.name_full:
        # Collection name
        coll_name = collection.name_full
        # Save the word before the underscore
        force_name = coll_name[0:coll_name.find('_')] + '_Forces'
        # Loop through the objects existed in the collection 
        for obj in collection.objects:
            # Cloth > Effector collection
            obj.modifiers['Cloth'].settings.effector_weights.collection = bpy.data.collections[force_name]
            # Cloth
            obj.modifiers['Cloth'].settings.quality = 6
            obj.modifiers['Cloth'].settings.air_damping = 5
            obj.modifiers['Cloth'].settings.bending_model = 'ANGULAR'
            # Cloth > Stiffness 
            obj.modifiers['Cloth'].settings.tension_stiffness = 2
            obj.modifiers['Cloth'].settings.compression_stiffness = 2
            obj.modifiers['Cloth'].settings.shear_stiffness = 2
            obj.modifiers['Cloth'].settings.bending_stiffness = 2
            # Cloth > Damping
            obj.modifiers['Cloth'].settings.tension_damping = 20
            obj.modifiers['Cloth'].settings.compression_damping = 20
            obj.modifiers['Cloth'].settings.shear_damping = 20
            obj.modifiers['Cloth'].settings.bending_damping = 0.5
            # Cloth > Internal Springs
            obj.modifiers['Cloth'].settings.use_internal_springs = True
            obj.modifiers['Cloth'].settings.internal_spring_max_length = 0
            obj.modifiers['Cloth'].settings.internal_spring_max_diversion = 0.785398
            obj.modifiers['Cloth'].settings.internal_spring_normal_check = True
            obj.modifiers['Cloth'].settings.internal_tension_stiffness = 1
            obj.modifiers['Cloth'].settings.internal_compression_stiffness = 1
            obj.modifiers['Cloth'].settings.internal_tension_stiffness_max = 5
            obj.modifiers['Cloth'].settings.internal_compression_stiffness_max = 5
            # Cloth > Pressure
            obj.modifiers['Cloth'].settings.uniform_pressure_force = 5
            obj.modifiers['Cloth'].settings.use_pressure_volume = True
            obj.modifiers['Cloth'].settings.target_volume = 1
            obj.modifiers['Cloth'].settings.pressure_factor = 1
            obj.modifiers['Cloth'].settings.fluid_density = 0
            # Cloth > Collisions
            obj.modifiers['Cloth'].collision_settings.collision_quality = 6
            obj.modifiers['Cloth'].collision_settings.use_collision = True
            obj.modifiers['Cloth'].collision_settings.use_self_collision = True
            obj.modifiers['Cloth'].collision_settings.self_friction = 5
            obj.modifiers['Cloth'].collision_settings.self_distance_min = 0.015
            obj.modifiers['Cloth'].collision_settings.self_impulse_clamp = 0
            # Collision
            obj.modifiers['Collision'].settings.damping = 0.579821
            obj.modifiers['Collision'].settings.thickness_outer = 0.02
            obj.modifiers['Collision'].settings.thickness_inner = 0.2
            obj.modifiers['Collision'].settings.cloth_friction = 5
            obj.modifiers['Collision'].settings.use_culling = True
            
handlers = goo.handler_class()
handlers.forces = [fA11, fA12, fA13, fA14, fA15, fA16, fA17, fA18, \
                    fA21, fA22, fA23, fA24, fA25, fA26, fA27, fA28, \
                    fA31, fA32, fA33, fA34, fA35, fA36, fA37, fA38, \
                    fA41, fA42, fA43, fA44, fA45, fA46, fA47, fA48, \
                    fA51, fA52, fA53, fA54, fA55, fA56, fA57, fA58]
bpy.app.handlers.frame_change_post.clear()
bpy.app.handlers.frame_change_post.append(handlers.adhesion_handler)
