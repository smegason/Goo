# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy
from importlib import reload

reload(goo)
goo.setup_world()

#================== Cell A Collection ==================
# Define cell A1
goo.make_cell("cell_A1", loc = (-0.95,0,0),  type = "cellsA")
# Define cell A2
goo.make_cell("cell_A2", loc = (0.95,0,0),  type = "cellsA")

# Define cell A1
goo.make_cell("cell_B3", loc = (-0.95,2.2,0), material = ('green', 0,.1,0), type = "cellsB")
# Define cell A2
goo.make_cell("cell_B4", loc = (0.95,2.2,0), material = ('green', 0,.1,0), type = "cellsB")


# Define cell A1
goo.make_cell("cell_C5", loc = (-0.95,4.4,0), material = ('red', .1,0,0), type = "cellsC")
# Define cell A2
goo.make_cell("cell_C6", loc = (0.95,4.4,0), material = ('red', .1,0,0), type = "cellsC")

#================== Force A Collection ==================

homoA = 750
homoB = 1500
homoC = 3000

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)
#goo.add_hetero_adhesion('cell_A1', 'cellsB', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_A2', -homoA)
#goo.add_hetero_adhesion('cell_A2', 'cellsB', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_B3', -homoB)
#goo.add_hetero_adhesion('cell_A3', 'cellsB', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_B4', -homoB)
#goo.add_hetero_adhesion('cell_A4', 'cellsB', -heteroAB)


# Define force A1
goo.add_homo_adhesion('cell_C5', -homoC)
#goo.add_hetero_adhesion('cell_B5', 'cellsA', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_C6', -homoC)
#goo.add_hetero_adhesion('cell_B6', 'cellsA', -heteroAB)


#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 500, # default, 250
                           filepath = "C:\\tmp\\sorting_test_rendering\\data5", 
                           adhesion = True, # default, True
                           data = False, # default, False
                           growth = True, 
                           division = False, 
                           motility = False
                           )