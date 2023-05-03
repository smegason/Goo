# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
import bpy

goo.setup_world()
  
#====================== Colors =========================
goo.add_material("green", 0, 0.1, 0)
goo.add_material("red", 0.1, 0, 0)

#================== Cell A Collection ==================
# Create a collection for cell A
goo.make_collection("A_Cells", type ='cell')
# Define cell A1
goo.make_cell("cell_A1", loc = (0,0,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A2", loc = (2,-2,0), collection = "A_Cells")
# Define cell A3
goo.make_cell("cell_A3", loc = (2,2,-0.5),  collection = "A_Cells")
# Define cell A4
goo.make_cell("cell_A4", loc = (0,0,2), collection = "A_Cells")
# Define cell A5
goo.make_cell("cell_A5", loc = (2,0,0.5), collection = "A_Cells")
# Define cell A6
goo.make_cell("cell_A6", loc = (0,2,0), collection = "A_Cells")
# Define cell A7
goo.make_cell("cell_A7", loc = (0,2,2), collection = "A_Cells")
# Define cell A8
goo.make_cell("cell_A8", loc = (-3,-2,1), collection = "A_Cells")
# Define cell A9
goo.make_cell("cell_A9", loc = (0,-2,-2), collection = "A_Cells")
# Define cell A10
goo.make_cell("cell_A10", loc = (-2,-2,-1), collection = "A_Cells")
# Define cell A11
goo.make_cell("cell_A11", loc = (-2,2,0.5), collection = "A_Cells")
# Define cell A12
goo.make_cell("cell_A12", loc = (0,-2,0), collection = "A_Cells")
# Define cell A13
goo.make_cell("cell_A13", loc = (2,4,1), collection = "A_Cells")
# Define cell A14
goo.make_cell("cell_A14", loc = (0,0,4), collection = "A_Cells")
# Define cell A15
goo.make_cell("cell_A15", loc = (-2,0,3), collection = "A_Cells")
# Define cell A16
goo.make_cell("cell_A16", loc = (-2,0,-0.5), collection = "A_Cells")



#================== Force A Collection ==================
# Force parameters; strength and decay power
force_strength = -3000
force_falloff = 0

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
# Define force A9
goo.make_force("force_A9", "cell_A9", force_strength, force_falloff, "A_Forces")
# Define force A10
goo.make_force("force_A10", "cell_A10", force_strength, force_falloff, "A_Forces")
# Define force A11
goo.make_force("force_A11", "cell_A11", force_strength, force_falloff, "A_Forces")
# Define force A12
goo.make_force("force_A12", "cell_A12", force_strength, force_falloff, "A_Forces")
# Define force A13
goo.make_force("force_A13", "cell_A13", force_strength, force_falloff, "A_Forces")
# Define force A14
goo.make_force("force_A14", "cell_A14", force_strength, force_falloff, "A_Forces")
# Define force A15
goo.make_force("force_A15", "cell_A15", force_strength, force_falloff, "A_Forces")
# Define force A16
goo.make_force("force_A16", "cell_A16", force_strength, force_falloff, "A_Forces")

#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 500, # default, 250
                           filepath = "C:\\Users\\anr9744\\Projects\\Goo\\data\\16clump", 
                           adhesion = True, # default, True
                           data = False, # default, False
                           growth = False
                           )
