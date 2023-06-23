# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy
from importlib import reload

reload(goo)
goo.setup_world(seed = 1)

#================== Cell A Collection ==================
# Define cell A1
goo.make_cell("cell_A1", loc = (-1,0,0), type = "cellsA")
# Define cell A2
goo.make_cell("cell_A2", loc = (-3,0,0), type = "cellsA")
# Define cell A1
goo.make_cell("cell_A3", loc = (-1,2,0), type = "cellsA")
# Define cell A2
goo.make_cell("cell_A4", loc = (-3,2,0), type = "cellsA")



# Define cell A1
goo.make_cell("cell_B5", loc = (1,0,0), material = ('red', 0.1,0,0), type = "cellsB")
# Define cell A2
goo.make_cell("cell_B6", loc = (3,0,0), material = ('red', 0.1,0,0), type = "cellsB")
# Define cell A1
goo.make_cell("cell_B7", loc = (1,2,0), material = ('red', 0.1,0,0), type = "cellsB")
# Define cell A2
goo.make_cell("cell_B8", loc = (3,2,0), material = ('red', 0.1,0,0), type = "cellsB")


#================== Force A Collection ==================

homoA = 1000
homoB = 1000
heteroAB = 3500

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)
goo.add_hetero_adhesion('cell_A1', 'cellsB', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_A2', -homoA)
goo.add_hetero_adhesion('cell_A2', 'cellsB', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_A3', -homoA)
goo.add_hetero_adhesion('cell_A3', 'cellsB', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_A4', -homoA)
goo.add_hetero_adhesion('cell_A4', 'cellsB', -heteroAB)

# Define force A1
goo.add_homo_adhesion('cell_B5', -homoB)
goo.add_hetero_adhesion('cell_B5', 'cellsA', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_B6', -homoB)
goo.add_hetero_adhesion('cell_B6', 'cellsA', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_B7', -homoB)
goo.add_hetero_adhesion('cell_B7', 'cellsA', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_B8', -homoB)
goo.add_hetero_adhesion('cell_B8', 'cellsA', -heteroAB)


#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 300, # default, 250
                           filepath = "C:\\tmp\\adhesion_paper\\hetero_homo\\hetero", 
                           adhesion = True, # default, True
                           data = True, # default, False
                           growth = True, 
                           division = False, 
                           motility = False
                           )