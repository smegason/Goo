# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy
from importlib import reload

reload(goo)
goo.setup_world(seed = 1)

#================== Cell A Collection ==================
# Define cell A1
goo.make_cell("cell_A1", loc = (0,10,0), type = "cellsA")
# Define cell A2
goo.make_cell("cell_A2", loc = (10,10,0), type = "cellsA")
# Define cell A1
goo.make_cell("cell_A3", loc = (0,0,0), type = "cellsA")
# Define cell A2
goo.make_cell("cell_A4", loc = (10,0,0), type = "cellsA")

# Define cell A1
#goo.make_cell("cell_B3", loc = (3,2,0), material = ('red', 0.1,0,0), type = "cellsB")
# Define cell A2
#goo.make_cell("cell_B4", loc = (5,2,0), material = ('red', 0.1,0,0), type = "cellsB")
# Define cell A1
#goo.make_cell("cell_B5", loc = (3,2,-2.2), material = ('red', 0.1,0,0), type = "cellsB")
# Define cell A2
#goo.make_cell("cell_B6", loc = (5,2,-2.2), material = ('red', 0.1,0,0), type = "cellsB")


#================== Force A Collection ==================

# Define force A1
goo.add_homo_adhesion("cell_A1", -2000)
goo.add_motion('cell_A1', -0)
# Define force A2
goo.add_homo_adhesion("cell_A2", -2000)
goo.add_motion('cell_A2', -500)
# Define force A1
goo.add_homo_adhesion("cell_A3", -2000)
goo.add_motion('cell_A3', -1000)
# Define force A2
goo.add_homo_adhesion("cell_A4", -2000)
goo.add_motion('cell_A4', -2000)


# Define force A1
'''goo.make_force("force_B3", "cell_B3", 'cellsB', -1000, 0)
goo.add_motion('cell_B3', -500)
# Define force A2
goo.make_force("force_B4", "cell_B4", 'cellsB', -1000, 0)
goo.add_motion('cell_B4', -500)
# Define force A1
goo.make_force("force_B5", "cell_B5", 'cellsB', -1000, 0)
goo.add_motion('cell_B5', -500)
# Define force A2
goo.make_force("force_B6", "cell_B6", 'cellsB', -1000, 0)
goo.add_motion('cell_B6', -500)'''

# Add random motion for cell type A
#goo.add_motion('cellsB', -1000)

#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 50, # default, 250
                           filepath = "C:\\tmp\\differential_motion\\data2", 
                           adhesion = True, # default, True
                           data = True, # default, False
                           growth = True, 
                           division = False, 
                           motility = True
                           )