# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy
from importlib import reload

reload(goo)
goo.setup_world()

#================== Cell A Collection ==================
# Define cell A1
goo.make_cell("cell_A1", loc = (3,0,0), type = "cellsA")
# Define cell A2
goo.make_cell("cell_A2", loc = (5,0,0), type = "cellsA")
# Define cell A1
goo.make_cell("cell_A3", loc = (3,0,-2.2), type = "cellsA")
# Define cell A2
goo.make_cell("cell_A4", loc = (5,0,-2.2), type = "cellsA")

# Define cell A1
goo.make_cell("cell_B3", loc = (3,2,0), material = ('red', 0.1,0,0), type = "cellsB")
# Define cell A2
goo.make_cell("cell_B4", loc = (5,2,0), material = ('red', 0.1,0,0), type = "cellsB")
# Define cell A1
goo.make_cell("cell_B5", loc = (3,2,-2.2), material = ('red', 0.1,0,0), type = "cellsB")
# Define cell A2
goo.make_cell("cell_B6", loc = (5,2,-2.2), material = ('red', 0.1,0,0), type = "cellsB")


#================== Force A Collection ==================

# Define force A1
goo.make_force("force_A1", "cell_A1", 'cellsA', -2000, 0)
goo.add_motion('cell_A1', -500)
# Define force A2
goo.make_force("force_A2", "cell_A2", 'cellsA', -2000, 0)
# Add random motion for cell type A
goo.add_motion('cell_A2', -500)
# Define force A1
goo.make_force("force_A3", "cell_A3", 'cellsA', -2000, 0)
goo.add_motion('cell_A3', -500)
# Define force A2
goo.make_force("force_A4", "cell_A4", 'cellsA', -2000, 0)
# Add random motion for cell type A
goo.add_motion('cell_A4', -500)


# Define force A1
goo.make_force("force_B3", "cell_B3", 'cellsB', -1000, 0)
goo.add_motion('cell_B3', -500)
# Define force A2
goo.make_force("force_B4", "cell_B4", 'cellsB', -1000, 0)
goo.add_motion('cell_B4', -500)
# Define force A1
goo.make_force("force_B5", "cell_B5", 'cellsB', -1000, 0)
goo.add_motion('cell_B5', -500)
# Define force A2
goo.make_force("force_B6", "cell_B6", 'cellsB', -1000, 0)
goo.add_motion('cell_B6', -500)

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