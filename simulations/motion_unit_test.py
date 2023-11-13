
from goo import goo
import bpy
from importlib import reload

reload(goo)
goo.setup_world(seed = 1)

#================== Cell A Collection ==================
# Define cell A1
goo.make_cell("cell_A1", loc = (0,10,0), type = "cellsA")


#================== Force A Collection ==================

# Define force A1
goo.add_homo_adhesion("cell_A1", -2000)
goo.add_motion('cell_A1', -1000, distribution = 'uniform', size = 0.1)


#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 50, # default, 250
                           filepath = "C:\\tmp\\differential_motion\\data2", 
                           adhesion = True, # default, True
                           data = True, # default, False
                           growth = True, 
                           division = False, 
                           motility = False
                           )