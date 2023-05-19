# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy
from importlib import reload

reload(goo)
goo.setup_world()

  
#================== Cell A Collection ==================
# Create a collection for cell A
goo.make_collection("A_Cells", type = 'cell')
# Define cell A1
#goo.make_cell("cell_A1", loc = (0,-5,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A2", loc = (0,-3,0), collection = "A_Cells")
# Define cell A3
goo.make_cell("cell_A3", loc = (0,-1,0), collection = "A_Cells")
# Define cell A4
goo.make_cell("cell_A4", loc = (0,1,0), collection = "A_Cells")
# Define cell A3
#goo.make_cell("cell_A5", loc = (0,3,0), collection = "A_Cells")
# Define cell A4
#goo.make_cell("cell_A6", loc = (0,5,0), collection = "A_Cells")

# Define cell A1
#goo.make_cell("cell_A7", loc = (-2,-5,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A8", loc = (-2,-3,0), collection = "A_Cells")
# Define cell A3
goo.make_cell("cell_A9", loc = (-2,-1,0), collection = "A_Cells")
# Define cell A4
goo.make_cell("cell_A10", loc = (-2,1,0), collection = "A_Cells")
# Define cell A3
#goo.make_cell("cell_A11", loc = (-2,3,0), collection = "A_Cells")
# Define cell A4
#goo.make_cell("cell_A12", loc = (-2,5,0), collection = "A_Cells")

# Define cell A1
#goo.make_cell("cell_A13", loc = (-4,-5,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A14", loc = (-4,-3,0), collection = "A_Cells")
# Define cell A1
goo.make_cell("cell_A15", loc = (-4,-1,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A16", loc = (-4,1,0), collection = "A_Cells")
# Define cell A1
#goo.make_cell("cell_A17", loc = (-4,3,0), collection = "A_Cells")
# Define cell A2
#goo.make_cell("cell_A18", loc = (-4,5,0), collection = "A_Cells")

'''# Define cell A1
goo.make_cell("cell_A19", loc = (4,-5,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A20", loc = (4,-3,0), collection = "A_Cells")
# Define cell A1
goo.make_cell("cell_A21", loc = (4,-1,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A22", loc = (4,1,0), collection = "A_Cells")
# Define cell A1
goo.make_cell("cell_A23", loc = (4,3,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A24", loc = (4,5,0), collection = "A_Cells")

# Define cell A1
goo.make_cell("cell_A25", loc = (2,-5,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A26", loc = (2,-3,0), collection = "A_Cells")
# Define cell A1
goo.make_cell("cell_A27", loc = (2,-1,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A28", loc = (2,1,0), collection = "A_Cells")
# Define cell A1
goo.make_cell("cell_A29", loc = (2,3,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A30", loc = (2,5,0), collection = "A_Cells")

# Define cell A1
goo.make_cell("cell_A31", loc = (6,-5,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A32", loc = (6,-3,0), collection = "A_Cells")
# Define cell A1
goo.make_cell("cell_A33", loc = (6,-1,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A34", loc = (6,1,0), collection = "A_Cells")
# Define cell A1
goo.make_cell("cell_A35", loc = (6,3,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A36", loc = (6,5,0), collection = "A_Cells")'''


#================== Force A Collection ==================

strength = -2000
falloff = 0

# Create a collection for force A
goo.make_collection("A_Forces", type = 'force')
# Define force A1
#goo.make_force("force_A1", "cell_A1", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A2", "cell_A2", strength, falloff, "A_Forces")
# Define force A3
goo.make_force("force_A3", "cell_A3", strength, falloff, "A_Forces")
# Define force A4
goo.make_force("force_A4", "cell_A4", strength, falloff, "A_Forces")
# Define force A5
#goo.make_force("force_A5", "cell_A5", strength, falloff, "A_Forces")
# Define force A6
#goo.make_force("force_A6", "cell_A6", strength, falloff, "A_Forces")

# Define force A1
#goo.make_force("force_A7", "cell_A7", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A8", "cell_A8", strength, falloff, "A_Forces")
# Define force A3
goo.make_force("force_A9", "cell_A9", strength, falloff, "A_Forces")
# Define force A4
goo.make_force("force_A10", "cell_A10", strength, falloff, "A_Forces")
# Define force A5
#goo.make_force("force_A11", "cell_A11", strength, falloff, "A_Forces")
# Define force A6
#goo.make_force("force_A12", "cell_A12", strength, falloff, "A_Forces")

# Define force A1
#goo.make_force("force_A13", "cell_A13", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A14", "cell_A14", strength, falloff, "A_Forces")
# Define force A1
goo.make_force("force_A15", "cell_A15", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A16", "cell_A16", strength, falloff, "A_Forces")
# Define force A1
#goo.make_force("force_A17", "cell_A17", strength, falloff, "A_Forces")
# Define force A2
#goo.make_force("force_A18", "cell_A18", strength, falloff, "A_Forces")

'''# Define force A1
goo.make_force("force_A19", "cell_A19", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A20", "cell_A20", strength, falloff, "A_Forces")
# Define force A1
goo.make_force("force_A21", "cell_A21", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A22", "cell_A22", strength, falloff, "A_Forces")
# Define force A1
goo.make_force("force_A23", "cell_A23", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A24", "cell_A24", strength, falloff, "A_Forces")

# Define force A1
goo.make_force("force_A25", "cell_A25", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A26", "cell_A26", strength, falloff, "A_Forces")
# Define force A1
goo.make_force("force_A27", "cell_A27", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A28", "cell_A28", strength, falloff, "A_Forces")
# Define force A1
goo.make_force("force_A29", "cell_A29", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A30", "cell_A30", strength, falloff, "A_Forces")

# Define force A1
goo.make_force("force_A31", "cell_A31", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A32", "cell_A32", strength, falloff, "A_Forces")
# Define force A1
goo.make_force("force_A33", "cell_A33", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A34", "cell_A34", strength, falloff, "A_Forces")
# Define force A1
goo.make_force("force_A35", "cell_A35", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A36", "cell_A36", strength, falloff, "A_Forces")
'''
#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 500, # default, 250
                           filepath = "C:\\Users\\anr9744\\Projects\\Goo\\data\\16clump", 
                           motion_strength = -500,
                           adhesion = True, # default, True
                           data = False, # default, False
                           growth = True, 
                           motility = True
                           )
