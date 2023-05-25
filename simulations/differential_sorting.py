# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
from importlib import reload
import bpy
reload(goo)
goo.setup_world()

#================== Cell A Collection ==================

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

# Define cell A1
goo.make_cell("cell_A7", loc = (0,0,2), type = "A_Cells")
# Define cell A2
goo.make_cell("cell_A8", loc = (4,0,2), type = "A_Cells")
# Define cell A3
goo.make_cell("cell_A9", loc = (2,2,2), type = "A_Cells")
# Define cell A3
goo.make_cell("cell_A11", loc = (2,-2,2), type = "A_Cells")

# Define cell A4
goo.make_cell("cell_A10", loc = (0,0,-2), type = "A_Cells")
# Define cell A1
goo.make_cell("cell_A1", loc = (4,0,-2), type = "A_Cells")
# Define cell A4
goo.make_cell("cell_A12", loc = (2,2,-2), type = "A_Cells")
# Define cell A4
goo.make_cell("cell_A13", loc = (2,-2,-2), type = "A_Cells")



# Define cell A1
goo.make_cell("cell_B13", loc = (0,0,0), material = ('red', 0.2, 0, 0), type = "B_Cells")
# Define cell A2
goo.make_cell("cell_B14", loc = (4,0,0), material = ('red', 0.2, 0, 0), type = "B_Cells")
# Define cell A1
goo.make_cell("cell_B15", loc = (2,2,0), material = ('red', 0.2, 0, 0), type = "B_Cells")
# Define cell A1
goo.make_cell("cell_B17", loc = (2,-2,0), material = ('red', 0.2, 0, 0), type = "B_Cells")

# Define cell A2
goo.make_cell("cell_B18", loc = (2,0,-2), material = ('red', 0.2, 0, 0), type = "B_Cells")
# Define cell A2
goo.make_cell("cell_B16", loc = (0,2,-2), material = ('red', 0.2, 0, 0), type = "B_Cells")
# Define cell A1
goo.make_cell("cell_B19", loc = (4,2,-2), material = ('red', 0.2, 0, 0), type = "B_Cells")
# Define cell A1
goo.make_cell("cell_B25", loc = (0,-2,-2), material = ('red', 0.2, 0, 0), type = "B_Cells")
# Define cell A1
goo.make_cell("cell_B26", loc = (4,-2,-2), material = ('red', 0.2, 0, 0), type = "B_Cells")

# Define cell A2
goo.make_cell("cell_B20", loc = (2,0,2), material = ('red', 0.2, 0, 0), type = "B_Cells")
# Define cell A1
goo.make_cell("cell_B21", loc = (0,2,2), material = ('red', 0.2, 0, 0), type = "B_Cells")
# Define cell A2
goo.make_cell("cell_B22", loc = (4,2,2), material = ('red', 0.2, 0, 0), type = "B_Cells")
# Define cell A1
goo.make_cell("cell_B23", loc = (0,-2,2), material = ('red', 0.2, 0, 0), type = "B_Cells")
# Define cell A2
goo.make_cell("cell_B24", loc = (4,-2,2), material = ('red', 0.2, 0, 0), type = "B_Cells")


#================== Force A Collection ==================

strength = -2000
falloff = 0

# Define force A1
goo.make_force("force_A1", "cell_A1", "A_Cells", strength, falloff)
goo.add_motion('cell_A1', -500)
# Define force A2
goo.make_force("force_A2", "cell_A2", "A_Cells", strength, falloff)
goo.add_motion('cell_A2', -500)
# Define force A3
goo.make_force("force_A3", "cell_A3", "A_Cells", strength, falloff)
goo.add_motion('cell_A3', -500)
# Define force A4
goo.make_force("force_A4", "cell_A4", "A_Cells", strength, falloff)
goo.add_motion('cell_A4', -500)
# Define force A5
goo.make_force("force_A5", "cell_A5", "A_Cells", strength, falloff)
goo.add_motion('cell_A5', -500)
# Define force A6
goo.make_force("force_A6", "cell_A6", "A_Cells", strength, falloff)
goo.add_motion('cell_A6', -500)

# Define force A1
goo.make_force("force_A7", "cell_A7", "A_Cells", strength, falloff)
goo.add_motion('cell_A7', -500)
# Define force A2
goo.make_force("force_A8", "cell_A8", "A_Cells", strength, falloff)
goo.add_motion('cell_A8', -500)
# Define force A3
goo.make_force("force_A9", "cell_A9", "A_Cells", strength, falloff)
goo.add_motion('cell_A9', -500)
# Define force A4
goo.make_force("force_A10", "cell_A10", "A_Cells", strength, falloff)
goo.add_motion('cell_A10', -500)
# Define force A5
goo.make_force("force_A11", "cell_A11", "A_Cells", strength, falloff)
goo.add_motion('cell_A11', -500)
# Define force A6
goo.make_force("force_A12", "cell_A12", "A_Cells", strength, falloff)
goo.add_motion('cell_A12', -500)
# Define force A6
goo.make_force("force_A13", "cell_A13", "A_Cells", strength, falloff)
goo.add_motion('cell_A13', -500)



# Define force A1
goo.make_force("force_B13", "cell_B13", "B_Cells", strength, falloff)
goo.add_motion('cell_B13', -1000)
# Define force A2
goo.make_force("force_B14", "cell_B14", "B_Cells", strength, falloff)
goo.add_motion('cell_B14', -1000)
# Define force A1
goo.make_force("force_B15", "cell_B15", "B_Cells", strength, falloff)
goo.add_motion('cell_B15', -1000)
# Define force A2
goo.make_force("force_B16", "cell_B16", "B_Cells", strength, falloff)
goo.add_motion('cell_B16', -1000)
# Define force A1
goo.make_force("force_B17", "cell_B17", "B_Cells", strength, falloff)
goo.add_motion('cell_B17', -1000)
# Define force A2
goo.make_force("force_B18", "cell_B18", "B_Cells", strength, falloff)
goo.add_motion('cell_B18', -1000)

# Define force A1
goo.make_force("force_B19", "cell_B19", "B_Cells", strength, falloff)
goo.add_motion('cell_B19', -1000)
# Define force A2
goo.make_force("force_B20", "cell_B20", "B_Cells", strength, falloff)
goo.add_motion('cell_B20', -1000)
# Define force A1
goo.make_force("force_B21", "cell_B21", "B_Cells", strength, falloff)
goo.add_motion('cell_B21', -1000)
# Define force A2
goo.make_force("force_B22", "cell_B22", "B_Cells", strength, falloff)
goo.add_motion('cell_B22', -1000)
# Define force A1
goo.make_force("force_B23", "cell_B23", "B_Cells", strength, falloff)
goo.add_motion('cell_B23', -1000)
# Define force A2
goo.make_force("force_B24", "cell_B24", "B_Cells", strength, falloff)
goo.add_motion('cell_B24', -1000)
# Define force A2
goo.make_force("force_B25", "cell_B25", "B_Cells", strength, falloff)
goo.add_motion('cell_B25', -1000)
# Define force A2
goo.make_force("force_B26", "cell_B26", "B_Cells", strength, falloff)
goo.add_motion('cell_B26', -1000)




#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, # default, 1
                           end = 1000, # default, 250
                           filepath = "C:\\Users\\anr9744\\Projects\\Goo\\data\\cell_sorting", 
                           adhesion = True, # default, True
                           data = False, # default, False
                           growth = True,
                           motility = True,
                           seed = 5
                           )