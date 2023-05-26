# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy
from importlib import reload

reload(goo)
goo.setup_world(seed = 1)

#================== Cell A Collection ==================
# Define cell A1
goo.make_cell("cell_A1", loc = (3,0,0), type = "cellsA", scale = (0.85, 0.85, 0.85))
# Define cell A2
goo.make_cell("cell_A2", loc = (5,0,0), type = "cellsA", scale = (0.75, 0.75, 0.75))

# Define cell A1
goo.make_cell("cell_B3", loc = (3,3,0), type = "cellsB")
# Define cell A2
goo.make_cell("cell_B4", loc = (5,3,0), type = "cellsB")



#================== Force A Collection ==================

# Define force A1
goo.make_force("force_A1", "cell_A1", 'cellsA', -2000, 0)
goo.add_motion('cell_A1', -500)

# Define force A2
goo.make_force("force_A2", "cell_A2", 'cellsA', -2000, 0)
# Add random motion for cell type A
goo.add_motion('cell_A2', -500)

# Define force A1
goo.make_force("force_B3", "cell_B3", 'cellsB', -1000, 0)
# Define force A2
goo.make_force("force_B4", "cell_B4", 'cellsB', -1000, 0)
# Add random motion for cell type A
#goo.add_motion('cellsB', -1000)

#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 500, # default, 250
                           filepath = "C:\\Users\\anr9744\\Projects\\Goo\\data\\cell_so1", 
                           adhesion = True, # default, True
                           data = True, # default, False
                           growth = True, 
                           division = False, 
                           motility = True
                           )