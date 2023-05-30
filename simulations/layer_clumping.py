# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
from importlib import reload
import bpy
reload(goo)
goo.setup_world(seed = 1)

#================== Cell A Collection ==================
# Create a collection for cell A
goo.make_collection("A_Cells", type = 'cell')
# Define cell A1
goo.make_cell("cell_A1", loc = (0,0,0), type = "A_Cells")
# Define cell A2
goo.make_cell("cell_A2", loc = (2,0,0), type = "A_Cells")
# Define cell A3
goo.make_cell("cell_A3", loc = (0,2,0), type = "A_Cells")
# Define cell A4
goo.make_cell("cell_A4", loc = (4,2,0), type = "A_Cells")
# Define cell A3
goo.make_cell("cell_A5", loc = (0,-2,0), type = "A_Cells")
# Define cell A4
goo.make_cell("cell_A6", loc = (4,-2,0), type = "A_Cells")
# Define cell A7
goo.make_cell("cell_A7", loc = (4,0,0), type = "A_Cells")
# Define cell A8
goo.make_cell("cell_A8", loc = (2,2,0), type = "A_Cells")
# Define cell A8
goo.make_cell("cell_A9", loc = (2,-2,0), type = "A_Cells")


#================== Force A Collection ==================
# Force parameters; strength and decay power
strength = -2000
falloff = 0
motion = -750

# Define force A1
goo.make_force("force_A1", "cell_A1", "A_Cells", strength, falloff)
goo.add_motion('cell_A1', motion)
# Define force A2
goo.make_force("force_A2", "cell_A2", "A_Cells", strength, falloff)
goo.add_motion('cell_A2', motion)
# Define force A3
goo.make_force("force_A3", "cell_A3", "A_Cells", strength, falloff)
goo.add_motion('cell_A3', motion)
# Define force A4
goo.make_force("force_A4", "cell_A4", "A_Cells", strength, falloff)
goo.add_motion('cell_A4', motion)
# Define force A5
goo.make_force("force_A5", "cell_A5", "A_Cells", strength, falloff)
goo.add_motion('cell_A5', motion)
# Define force A6
goo.make_force("force_A6", "cell_A6", "A_Cells", strength, falloff)
goo.add_motion('cell_A6', motion)

# Define force A1
goo.make_force("force_A7", "cell_A7", "A_Cells", strength, falloff)
goo.add_motion('cell_A7', motion)
# Define force A2
goo.make_force("force_A8", "cell_A8", "A_Cells", strength, falloff)
goo.add_motion('cell_A8', motion)
# Define force A3
goo.make_force("force_A9", "cell_A9", "A_Cells", strength, falloff)
goo.add_motion('cell_A9', motion)

#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 500, # default, 250
                           filepath = "C:\\tmp\\layer_clumping", 
                           adhesion = True, # default, True
                           data = False, # default, False
                           growth = True, 
                           motility = True
                           )