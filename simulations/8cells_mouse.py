# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy
goo.setup_world()


#================== Cell A Collection ==================
# Create a collection for cell A
goo.make_collection("A_Cells", type = 'cell')
# Define cell A1
goo.make_cell("cell_A1", loc = (-1.36,-1,0.2), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A2", loc = (-0.2,-1,-1.36), collection = "A_Cells")
# Define cell A3
goo.make_cell("cell_A3", loc = (1.36,-1,-0.2), collection = "A_Cells")
# Define cell A4
goo.make_cell("cell_A4", loc = (0.2,-1,1.36), collection = "A_Cells")
# Define cell A5
goo.make_cell("cell_A5", loc = (-0.73,1,1.1), collection = "A_Cells")
# Define cell A6
goo.make_cell("cell_A6", loc = (1.2,1,0.84), collection = "A_Cells")
# Define cell A7
goo.make_cell("cell_A7", loc = (0.92,1,-1.08), collection = "A_Cells")
# Define cell A8
goo.make_cell("cell_A8", loc = (-1,1,-0.82), collection = "A_Cells")


#================== Force A Collection ==================
# Force parameters; strength and decay power
force_strength = -1200
force_falloff = 1

# Create a collection for force A
goo.make_collection("A_Forces", type = 'force')
# Define force A1
goo.make_force("force_A1", "cell_A1", force_strength, force_falloff, "A_Forces")
# Define force A2
goo.make_force("force_A2", "cell_A2", force_strength, force_falloff, "A_Forces")
# Define force A3
goo.make_force("force_A3", "cell_A3", force_strength, force_falloff, "A_Forces")
# Define force A4
goo.make_force("force_A4", "cell_A4", force_strength, force_falloff, "A_Forces")
# Define force A5
goo.make_force("force_A5", "cell_A5", force_strength, force_falloff, "A_Forces")
# Define force A6
goo.make_force("force_A6", "cell_A6", force_strength, force_falloff, "A_Forces")
# Define force A7
goo.make_force("force_A7", "cell_A7", force_strength, force_falloff, "A_Forces")
# Define force A8
goo.make_force("force_A8", "cell_A8", force_strength, force_falloff, "A_Forces")


#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 500, # default, 250
                           filepath = "C:\\Users\\anr9744\\Projects\\Goo\\data\\TEST_doublet_simulation", 
                           adhesion = True, # default, True
                           data = False, # default, False
                           )