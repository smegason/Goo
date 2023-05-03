# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy
from importlib import reload

reload(goo)
goo.setup_world()

#================== Cell A Collection ==================
# Create a collection for cell A
goo.make_collection(name = "cellsA", type = 'cell')
# Define cell A1
goo.make_cell("cell_A1", loc = (3,0,0), collection = "cellsA", scale = (0.85, 0.85, 0.85))
# Define cell A2
goo.make_cell("cell_A2", loc = (5,0,0), collection = "cellsA", scale = (0.75, 0.75, 0.75))

# Create a collection for cell A
'''goo.make_collection(name = "cellsB", type = 'cell')
# Define cell A1
goo.make_cell("cell_A3", loc = (3,3,0), collection = "cellsB")
# Define cell A2
goo.make_cell("cell_A4", loc = (5,3,0), collection = "cellsB")
'''


#================== Force A Collection ==================

# Create a collection for force A
goo.make_collection(name = "forcesA", type = 'force')         
# Define force A1
goo.make_force("force_A1", "cell_A1", -2000, 0, "forcesA")
# Define force A2
goo.make_force("force_A2", "cell_A2", -2000, 0, "forcesA")

# Create a collection for force A
'''goo.make_collection(name = "forcesB", type = 'force')         
# Define force A1
goo.make_force("force_A3", "cell_A3", -1000, 1, "forcesB")
# Define force A2
goo.make_force("force_A4", "cell_A4", -1000, 1, "forcesB")'''

#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 500, # default, 250
                           filepath = "C:\\Users\\anr9744\\Projects\\Goo\\data\\cell_so1", 
                           adhesion = True, # default, True
                           data = True, # default, False
                           growth = True, 
                           division = False
                           )