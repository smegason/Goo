# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy
from importlib import reload

reload(goo)
goo.setup_world()

#================== Cell A Collection ==================
# Create a collection for cell A
goo.make_collection("cells", 'cell')
# Define cell A1
goo.make_cell("cell_A1", stiffness = 10, loc = (0,0,0), collection = "cells", scale = (1,1,1))
# Define cell A2
#goo.make_cell("cell_A2", loc = (5,0,0), collection = "A_Cells")


#================== Force A Collection ==================

# Create a collection for force A
goo.make_collection("forces", 'force')         
# Define force A1
goo.make_force("force_A1", "cell_A1", -1000, 1, "forces")
# Define force A2
#goo.make_force("force_A2", "cell_A2", -1000, 1, "A_Forces")

#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 500, # default, 250
                           filepath = "C:\\Users\\anr9744\\Projects\\Goo\\data\\cell_so", 
                           adhesion = True, # default, True
                           data = False, # default, False
                           growth = False, 
                           division = True
                           )