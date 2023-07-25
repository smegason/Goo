# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy
from importlib import reload

reload(goo)
goo.setup_world()

#================== Cell A Collection ==================
# Define cell A1
goo.make_cell("cell_A1", loc = (-0.34,-0.15,-0.86), type = "cellsA")
# Define cell A2
goo.make_cell("cell_A2", loc = (2.76,0,0.56), type = "cellsA")
# Define cell A1
goo.make_cell("cell_A3", loc = (-1,2.69,-0.15), type = "cellsA")
# Define cell A2
goo.make_cell("cell_A4", loc = (0.2,0.4,1.39), type = "cellsA")


# Define cell A1
goo.make_cell("cell_B5", loc = (0.69,2.88,0.93), type = "cellsA")
# Define cell A2
goo.make_cell("cell_B6", loc = (2,2.66,-0.57), type = "cellsA")
# Define cell A1
goo.make_cell("cell_B7", loc = (1.31,1.16,-1.96), type = "cellsA")
# Define cell A2
#goo.make_cell("cell_B8", loc = (2.08,1.5,1.87), type = "cellsA")


#================== Force A Collection ==================

homoA = 2500
homoB = 2500
#heteroAB = 500
motion = 1000

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)
#goo.add_hetero_adhesion('cell_A1', 'cellsB', -heteroAB)
goo.add_motion('cell_A1', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A2', -homoA)
#goo.add_hetero_adhesion('cell_A2', 'cellsB', -heteroAB)
goo.add_motion('cell_A2', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A3', -homoA)
#goo.add_hetero_adhesion('cell_A3', 'cellsB', -heteroAB)
goo.add_motion('cell_A3', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A4', -homoA)
#goo.add_hetero_adhesion('cell_A4', 'cellsB', -heteroAB)
goo.add_motion('cell_A4', -motion)


# Define force A1
goo.add_homo_adhesion('cell_B5', -homoB)
#goo.add_hetero_adhesion('cell_B5', 'cellsA', -heteroAB)
goo.add_motion('cell_B5', -motion)
# Define force A2
goo.add_homo_adhesion('cell_B6', -homoB)
#goo.add_hetero_adhesion('cell_B6', 'cellsA', -heteroAB)
goo.add_motion('cell_B6', -motion)
# Define force A2
goo.add_homo_adhesion('cell_B7', -homoB)
#goo.add_hetero_adhesion('cell_B7', 'cellsA', -heteroAB)
goo.add_motion('cell_B7', -motion)
# Define force A2
#goo.add_homo_adhesion('cell_B8', -homoB)
#goo.add_hetero_adhesion('cell_B8', 'cellsA', -heteroAB)
#goo.add_motion('cell_B8', -motion)

#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 500, # default, 250
                           filepath = "C:\\tmp\\sorting_test_rendering\\data5", 
                           adhesion = True, # default, True
                           data = False, # default, False
                           growth = True, 
                           division = False, 
                           motility = True
                           )