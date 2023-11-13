# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy
from importlib import reload
reload(goo)
goo.setup_world()

  
#================== Cell A Collection ==================
# Define cell A1
goo.make_cell("cell_A1", loc = (0,0,0), type = "A_Cells")

#================== Force A Collection ==================

homoA = 1000

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)

#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 500, # default, 250
                           filepath = "C:\\Users\\tmp\\", 
                           adhesion = True, # default, True
                           data = False, # default, False
                           growth = True, 
                           motility = False, 
                           division = True
                           )
