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

# Define cell A2
cA2 = goo.Cell("cell_A2", loc = (2,-2,0))
# Make a Blender mesh object for cell
goo.make_cell(cA2)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A3
cA3 = goo.Cell("cell_A3", loc = (4,-2,0))
# Make a Blender mesh object for cell
goo.make_cell(cA3)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A4
cA4 = goo.Cell("cell_A4", loc = (4,-4,2))
# Make a Blender mesh object for cell
goo.make_cell(cA4)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A5
cA5 = goo.Cell("cell_A5", loc = (2,2,2))
# Make a Blender mesh object for cell
goo.make_cell(cA5)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A6
cA6 = goo.Cell("cell_A6", loc = (1,1,2))
# Make a Blender mesh object for cell
goo.make_cell(cA6)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A7
cA7 = goo.Cell("cell_A7", loc = (4,2,4))
# Make a Blender mesh object for cell
goo.make_cell(cA7)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A8
cA8 = goo.Cell("cell_A8", loc = (-1,-1,-2))
# Make a Blender mesh object for cell
goo.make_cell(cA8)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A9
cA9 = goo.Cell("cell_A9", loc = (5,2,-2))
# Make a Blender mesh object for cell
goo.make_cell(cA9)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene 
# if you do not delete them the object outside colleciton will link to object inside the collection
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A10
cA10 = goo.Cell("cell_A10", loc = (-5,-2,-3))
# Make a Blender mesh object for cell
goo.make_cell(cA10)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A11
cA11 = goo.Cell("cell_A11", loc = (4,-5,0))
# Make a Blender mesh object for cell
goo.make_cell(cA11)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A12
cA12 = goo.Cell("cell_A12", loc = (-3,-4,1))
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

# Define cell A13
cA13 = goo.Cell("cell_A13", loc = (6,6,6))
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
cA14 = goo.Cell("cell_A14", loc = (8,1,6))
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
cA15 = goo.Cell("cell_A15", loc = (-7,-6,4))
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
cA16 = goo.Cell("cell_A16", loc = (-1,-4,-3))
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
cA17 = goo.Cell("cell_A17", loc = (-5,-4,-3))
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
'''
# Define cell A18
cA18 = goo.Cell("cell_A18", loc = (-3,4,6))
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
'''
# Define cell A19
cA19 = goo.Cell("cell_A19", loc = (-3,-4,-1))
# Make a Blender mesh object for cell
goo.make_cell(cA19)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A20
cA20 = goo.Cell("cell_A20", loc = (3,3,-6))
# Make a Blender mesh object for cell
goo.make_cell(cA20)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A21
cA21 = goo.Cell("cell_A21", loc = (4,1,4))
# Make a Blender mesh object for cell
goo.make_cell(cA21)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A22
cA22 = goo.Cell("cell_A22", loc = (-2,-6,4))
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
cA23 = goo.Cell("cell_A23", loc = (-1,-4,-5))
# Make a Blender mesh object for cell
goo.make_cell(cA23)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A24
cA24 = goo.Cell("cell_A24", loc = (4,-5,-4))
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
cA25 = goo.Cell("cell_A25", loc = (-6,-4,3))
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

# Define cell A26
cA26 = goo.Cell("cell_A26", loc = (2,-6,4))
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
cA27 = goo.Cell("cell_A27", loc = (6,-5,-5))
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
cA28 = goo.Cell("cell_A28", loc = (-8,-6,4))
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

# Define cell A29
cA29 = goo.Cell("cell_A29", loc = (1,-7,-5))
# Make a Blender mesh object for cell
goo.make_cell(cA29)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A30
cA30 = goo.Cell("cell_A30", loc = (-8,-6,-4))
# Make a Blender mesh object for cell
goo.make_cell(cA30)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A31
cA31 = goo.Cell("cell_A31", loc = (-1,-7,5))
# Make a Blender mesh object for cell
goo.make_cell(cA31)
# The created cell is the active object
obj = bpy.context.active_object
# Set the object color
obj.active_material = matg
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active cell to our specific collection 
bpy.data.collections['A_Cells'].objects.link(obj)

# Define cell A32
cA32 = goo.Cell("cell_A32", loc = (-1,4,-5))
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
cA33 = goo.Cell("cell_A33", loc = (6,4,-5))
# Make a Blender mesh object for cell
goo.make_cell(cA33)
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
# Define and link force A1
fA1 = goo.Force("force_A1", "cell_A1", -800, 1)
# Make force
goo.make_force(fA1)
# The created force is the active object
obj = bpy.context.active_object

# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A2
fA2 = goo.Force("force_A2", "cell_A2", -800, 1)
# Make force
goo.make_force(fA2)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A3
fA3 = goo.Force("force_A3", "cell_A3", -800, 1)
# Make force
goo.make_force(fA3)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A4
fA4 = goo.Force("force_A4", "cell_A4", -800, 1)
# Make force
goo.make_force(fA4)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A5
fA5 = goo.Force("force_A5", "cell_A5", -800, 1)
# Make force
goo.make_force(fA5)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A6
fA6 = goo.Force("force_A6", "cell_A6", -800, 1)
# Make force
goo.make_force(fA6)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A7
fA7 = goo.Force("force_A7", "cell_A7", -800, 1)
# Make force
goo.make_force(fA7)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A8
fA8 = goo.Force("force_A8", "cell_A8", -800, 1)
# Make force
goo.make_force(fA8)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A9
fA9 = goo.Force("force_A9", "cell_A9", -800, 1)
# Make force
goo.make_force(fA9)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A10
fA10 = goo.Force("force_A10", "cell_A10", -800, 1)
# Make force
goo.make_force(fA10)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A11
fA11 = goo.Force("force_A11", "cell_A11", -800, 1)
# Make force
goo.make_force(fA11)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A12
fA12 = goo.Force("force_A12", "cell_A12", -800, 1)
# Make force
goo.make_force(fA12)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A13
fA13 = goo.Force("force_A13", "cell_A13", -800, 1)
# Make force
goo.make_force(fA13)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A14
fA14 = goo.Force("force_A14", "cell_A14", -800, 1)
# Make force
goo.make_force(fA14)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A15
fA15 = goo.Force("force_A15", "cell_A15", -800, 1)
# Make force
goo.make_force(fA15)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A16
fA16 = goo.Force("force_A16", "cell_A16", -800, 1)
# Make force
goo.make_force(fA16)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A17
fA17 = goo.Force("force_A17", "cell_A17", -800, 1)
# Make force
goo.make_force(fA17)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)
'''
# Define and link force A18
fA18 = goo.Force("force_A18", "cell_A18", -800, 1)
# Make force
goo.make_force(fA18)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)
'''
# Define and link force A19
fA19 = goo.Force("force_A19", "cell_A19", -800, 1)
# Make force
goo.make_force(fA19)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A20
fA20 = goo.Force("force_A20", "cell_A20", -800, 1)
# Make force
goo.make_force(fA20)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A21
fA21 = goo.Force("force_A21", "cell_A21", -800, 1)
# Make force
goo.make_force(fA21)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A22
fA22 = goo.Force("force_A22", "cell_A22", -800, 1)
# Make force
goo.make_force(fA22)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A23
fA23 = goo.Force("force_A23", "cell_A23", -800, 1)
# Make force
goo.make_force(fA23)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A24
fA24 = goo.Force("force_A24", "cell_A24", -800, 1)
# Make force
goo.make_force(fA24)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A25
fA25 = goo.Force("force_A25", "cell_A25", -800, 1)
# Make force
goo.make_force(fA25)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A26
fA26 = goo.Force("force_A26", "cell_A26", -800, 1)
# Make force
goo.make_force(fA26)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A27
fA27 = goo.Force("force_A27", "cell_A27", -800, 1)
# Make force
goo.make_force(fA27)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A28
fA28 = goo.Force("force_A28", "cell_A28", -800, 1)
# Make force
goo.make_force(fA28)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A29
fA29 = goo.Force("force_A29", "cell_A29", -800, 1)
# Make force
goo.make_force(fA29)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A30
fA30 = goo.Force("force_A30", "cell_A30", -800, 1)
# Make force
goo.make_force(fA30)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A31
fA31 = goo.Force("force_A31", "cell_A31", -800, 1)
# Make force
goo.make_force(fA31)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A30
fA32 = goo.Force("force_A32", "cell_A32", -800, 1)
# Make force
goo.make_force(fA32)
# The created force is the active object
obj = bpy.context.active_object
# Remove object from all collections not used in a scene
bpy.ops.collection.objects_remove_all()
# Add the active force to our specific collection 
bpy.data.collections['A_Forces'].objects.link(obj)

# Define and link force A33
fA33 = goo.Force("force_A33", "cell_A33", -800, 1)
# Make force
goo.make_force(fA33)
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
            obj.modifiers['Cloth'].settings.tension_stiffness = 10
            obj.modifiers['Cloth'].settings.compression_stiffness = 10
            obj.modifiers['Cloth'].settings.shear_stiffness = 10
            obj.modifiers['Cloth'].settings.bending_stiffness = 10
            # Cloth > Damping
            obj.modifiers['Cloth'].settings.tension_damping = 1
            obj.modifiers['Cloth'].settings.compression_damping = 1
            obj.modifiers['Cloth'].settings.shear_damping = 1
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
handlers.forces = [fA1, fA2, fA3, fA4, fA5, fA6, fA7, fA8, fA9, fA10, fA11, fA12, fA13, \
                    fA14, fA15, fA16, fA17, fA19, fA20, fA21, fA22, fA23, fA24, fA25, fA26, \
                    fA27, fA28, fA29, fA30, fA31, fA32, fA33]
bpy.app.handlers.frame_change_post.clear()
bpy.app.handlers.frame_change_post.append(handlers.adhesion_handler)
