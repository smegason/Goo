# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
from importlib import reload
import bpy
reload(goo)
goo.setup_world()

#================== Cell A Collection ==================
# Create a collection for cell A
goo.make_collection("A_Cells", type = 'cell')
# Define cell A1
goo.make_cell("cell_A1", loc = (0.1,0,0), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A2", loc = (2,0,0), collection = "A_Cells")
# Define cell A3
goo.make_cell("cell_A3", loc = (0,1,1.83), collection = "A_Cells")
# Define cell A4
goo.make_cell("cell_A4", loc = (0.74,2.85,-0.57), collection = "A_Cells")
# Define cell A3
goo.make_cell("cell_A5", loc = (0,-2,0), collection = "A_Cells")
# Define cell A4
goo.make_cell("cell_A6", loc = (-2,0,0), collection = "A_Cells")
# Define cell A1
goo.make_cell("cell_A7", loc = (3.39,0.25,0.76), collection = "A_Cells")
# Define cell A2
goo.make_cell("cell_A8", loc = (2,-2.03,0), collection = "A_Cells")
# Define cell A3
goo.make_cell("cell_A9", loc = (0.098,-2.59,1.79), collection = "A_Cells")
# Define cell A4
goo.make_cell("cell_A10", loc = (1.67,0.73,2.9), collection = "A_Cells")
# Define cell A3
goo.make_cell("cell_A11", loc = (0.86,-1.22,-1.92), collection = "A_Cells")
# Define cell A4
goo.make_cell("cell_A12", loc = (-1.15,1,-1.93), collection = "A_Cells")

# Create a collection for cell A
goo.make_collection("B_Cells", type = 'cell')
# Define cell A1
goo.make_cell("cell_A13", loc = (2,2,1.4), material = ('red', 0.2, 0, 0), collection = "B_Cells")
# Define cell A2
goo.make_cell("cell_A14", loc = (-0.98,1.95,0), material = ('red', 0.2, 0, 0), collection = "B_Cells")
# Define cell A1
goo.make_cell("cell_A15", loc = (0.97,-0.59,1.63), material = ('red', 0.2, 0, 0), collection = "B_Cells")
# Define cell A2
goo.make_cell("cell_A16", loc = (2.66,2,-0.57), material = ('red', 0.2, 0, 0), collection = "B_Cells")
# Define cell A1
goo.make_cell("cell_A17", loc = (0.92,1.16,-1.47), material = ('red', 0.2, 0, 0), collection = "B_Cells")
# Define cell A2
goo.make_cell("cell_A18", loc = (-1.16,-1,-1.22), material = ('red', 0.2, 0, 0), collection = "B_Cells")

# Define cell A1
goo.make_cell("cell_A19", loc = (-1.07,-0.76,1.58), material = ('red', 0.2, 0, 0), collection = "B_Cells")
# Define cell A2
goo.make_cell("cell_A20", loc = (0.1,2.9,1.38), material = ('red', 0.2, 0, 0), collection = "B_Cells")
# Define cell A1
goo.make_cell("cell_A21", loc = (1.9,-2.45,1.9), material = ('red', 0.2, 0, 0), collection = "B_Cells")
# Define cell A2
goo.make_cell("cell_A22", loc = (-0.05,-1.18,3.44), material = ('red', 0.2, 0, 0), collection = "B_Cells")
# Define cell A1
goo.make_cell("cell_A23", loc = (2.74,-0.34,-1.8), material = ('red', 0.2, 0, 0), collection = "B_Cells")
# Define cell A2
goo.make_cell("cell_A24", loc = (-2.1,0.9,1.81), material = ('red', 0.2, 0, 0), collection = "B_Cells")




#================== Force A Collection ==================

strength = -2000
falloff = 0

# Create a collection for force A
goo.make_collection("A_Forces", type = 'force')
# Define force A1
goo.make_force("force_A1", "cell_A1", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A2", "cell_A2", strength, falloff, "A_Forces")
# Define force A3
goo.make_force("force_A3", "cell_A3", strength, falloff, "A_Forces")
# Define force A4
goo.make_force("force_A4", "cell_A4", strength, falloff, "A_Forces")
# Define force A5
goo.make_force("force_A5", "cell_A5", strength, falloff, "A_Forces")
# Define force A6
goo.make_force("force_A6", "cell_A6", strength, falloff, "A_Forces")

# Define force A1
goo.make_force("force_A7", "cell_A7", strength, falloff, "A_Forces")
# Define force A2
goo.make_force("force_A8", "cell_A8", strength, falloff, "A_Forces")
# Define force A3
goo.make_force("force_A9", "cell_A9", strength, falloff, "A_Forces")
# Define force A4
goo.make_force("force_A10", "cell_A10", strength, falloff, "A_Forces")
# Define force A5
goo.make_force("force_A11", "cell_A11", strength, falloff, "A_Forces")
# Define force A6
goo.make_force("force_A12", "cell_A12", strength, falloff, "A_Forces")


# Create a collection for force A
goo.make_collection("B_Forces", type = 'force')
# Define force A1
goo.make_force("force_A13", "cell_A13", strength, falloff, "B_Forces")
# Define force A2
goo.make_force("force_A14", "cell_A14", strength, falloff, "B_Forces")
# Define force A1
goo.make_force("force_A15", "cell_A15", strength, falloff, "B_Forces")
# Define force A2
goo.make_force("force_A16", "cell_A16", strength, falloff, "B_Forces")
# Define force A1
goo.make_force("force_A17", "cell_A17", strength, falloff, "B_Forces")
# Define force A2
goo.make_force("force_A18", "cell_A18", strength, falloff, "B_Forces")

# Define force A1
goo.make_force("force_A19", "cell_A19", strength, falloff, "B_Forces")
# Define force A2
goo.make_force("force_A20", "cell_A20", strength, falloff, "B_Forces")
# Define force A1
goo.make_force("force_A21", "cell_A21", strength, falloff, "B_Forces")
# Define force A2
goo.make_force("force_A22", "cell_A22", strength, falloff, "B_Forces")
# Define force A1
goo.make_force("force_A23", "cell_A23", strength, falloff, "B_Forces")
# Define force A2
goo.make_force("force_A24", "cell_A24", strength, falloff, "B_Forces")


#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 500, # default, 250
                           filepath = "C:\\Users\\anr9744\\Projects\\Goo\\data\\cell_sorting", 
                           adhesion = True, # default, True
                           data = False, # default, False
                           growth = True, 
                           motility = True
                           )