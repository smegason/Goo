# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy
from importlib import reload
reload(goo)
goo.setup_world()

  
#================== Cell A Collection ==================
# Create a type for cell A
#goo.make_type("A_Cells", type = 'cell')
# Define cell A1
goo.make_cell("cell_A1", loc = (-3,1.5,1.5), type = "A_Cells")
# Define cell A2
goo.make_cell("cell_A2", loc = (-3,-1.5,1.5), type = "A_Cells")
# Define cell A3
goo.make_cell("cell_A3", loc = (0,1.5,1.5), type = "A_Cells")
# Define cell A4
goo.make_cell("cell_A4", loc = (0,-1.5,1.5), type = "A_Cells")
# Define cell A3
goo.make_cell("cell_A5", loc = (3,1.5,1.5), type = "A_Cells")
# Define cell A4
goo.make_cell("cell_A6", loc = (3,-1.5,1.5), type = "A_Cells")
# Define cell A1
goo.make_cell("cell_A7",  loc = (-3,1.5,-1.5), type = "A_Cells")
# Define cell A2
goo.make_cell("cell_A8", loc = (-3,-1.5,-1.5), type = "A_Cells")
# Define cell A3
goo.make_cell("cell_A9", loc = (0,1.5,-1.5), type = "A_Cells")
# Define cell A4
goo.make_cell("cell_A10", loc = (0,-1.5,-1.5), type = "A_Cells")
# Define cell A3
goo.make_cell("cell_A11", loc = (3,1.5,-1.5), type = "A_Cells")
# Define cell A4
goo.make_cell("cell_A12", loc = (3,-1.5,-1.5), type = "A_Cells")


#================== Force A Collection ==================

homoA = 2500
motion = 500

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)
goo.add_motion('cell_A1', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A2', -homoA)
goo.add_motion('cell_A2', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A3', -homoA)
goo.add_motion('cell_A3', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A4', -homoA)
goo.add_motion('cell_A4', -motion)
# Define force A1
goo.add_homo_adhesion('cell_A5', -homoA)
goo.add_motion('cell_A5', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A6', -homoA)
goo.add_motion('cell_A6', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A7', -homoA)
goo.add_motion('cell_A7', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A8', -homoA)
goo.add_motion('cell_A8', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A9', -homoA)
goo.add_motion('cell_A9', -motion)
# Define force A1
goo.add_homo_adhesion('cell_A10', -homoA)
goo.add_motion('cell_A10', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A11', -homoA)
goo.add_motion('cell_A11', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A12', -homoA)
goo.add_motion('cell_A12', -motion)

#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 500, # default, 250
                           filepath = "C:\\Users\\tmp\\", 
                           adhesion = True, # default, True
                           data = True, # default, False
                           growth = True, 
                           motility = True
                           )
