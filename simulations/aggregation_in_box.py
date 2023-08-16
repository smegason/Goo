# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy
from importlib import reload

reload(goo)
goo.setup_world()

#================== Cell A Collection ==================
# Define cell A1
goo.make_cell("cell_A1", loc = (-3,1.5,1.5), type = "cellsA", stiffness = 1)
# Define cell A2
goo.make_cell("cell_A2", loc = (-3,-1.5,1.5), type = "cellsA", stiffness = 1)
# Define cell A3
goo.make_cell("cell_A3", loc = (0,1.5,1.5), type = "cellsA", stiffness = 1)
# Define cell A4
goo.make_cell("cell_A4", loc = (0,-1.5,1.5), type = "cellsA", stiffness = 1)
# Define cell A3
goo.make_cell("cell_A5", loc = (0,1.5,-1.5), type = "cellsA", stiffness = 1)
# Define cell A4
goo.make_cell("cell_A6", loc = (0,-1.5,-1.5), type = "cellsA", stiffness = 1)
# Define cell A1
goo.make_cell("cell_A7",  loc = (-3,1.5,-1.5), type = "cellsA", stiffness = 1)
# Define cell A2
goo.make_cell("cell_A8", loc = (-3,-1.5,-1.5), type = "cellsA", stiffness = 1)

#================== Force A Collection ==================

homoA = 2000
#heteroAB = 500
motion = 500
size = 2

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)
#goo.add_hetero_adhesion('cell_A1', 'cellsB', -heteroAB)
goo.add_motion('cell_A1', -motion, size = size)
# Define force A2
goo.add_homo_adhesion('cell_A2', -homoA)
#goo.add_hetero_adhesion('cell_A2', 'cellsB', -heteroAB)
goo.add_motion('cell_A2', -motion, size = size)
# Define force A2
goo.add_homo_adhesion('cell_A3', -homoA)
#goo.add_hetero_adhesion('cell_A3', 'cellsB', -heteroAB)
goo.add_motion('cell_A3', -motion, size = size)
# Define force A2
goo.add_homo_adhesion('cell_A4', -homoA)
#goo.add_hetero_adhesion('cell_A4', 'cellsB', -heteroAB)
goo.add_motion('cell_A4', -motion, size = size)


# Define force A1
goo.add_homo_adhesion('cell_A5', -homoA)
#goo.add_hetero_adhesion('cell_B5', 'cellsA', -heteroAB)
goo.add_motion('cell_A5', -motion, size = size)
# Define force A2
goo.add_homo_adhesion('cell_A6', -homoA)
#goo.add_hetero_adhesion('cell_B6', 'cellsA', -heteroAB)
goo.add_motion('cell_A6', -motion, size = size)
# Define force A2
goo.add_homo_adhesion('cell_A7', -homoA)
#goo.add_hetero_adhesion('cell_B7', 'cellsA', -heteroAB)
goo.add_motion('cell_A7', -motion, size = size)
# Define force A2
goo.add_homo_adhesion('cell_A8', -homoA)
#goo.add_hetero_adhesion('cell_B8', 'cellsA', -heteroAB)
goo.add_motion('cell_A8', -motion, size = size)


#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 5000, # default, 250
                           filepath = "C:\\tmp\\aggregation\\data4",
                           adhesion = True, # default, True
                           data = True, # default, False
                           growth = True, 
                           division = False, 
                           motility = True
                           )
                           
                           
bpy.ops.mesh.primitive_cube_add(enter_editmode=False, align='WORLD', location=(-1.5, 0, 0), scale=(1, 1, 1))
bpy.context.object.scale[0] = -4.5
bpy.context.object.scale[1] = -4.5
bpy.context.object.scale[2] = -4.5
bpy.ops.object.modifier_add(type='COLLISION')
