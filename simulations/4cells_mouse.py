
from goo import goo
import bpy
from importlib import reload

reload(goo)
goo.setup_world()

#================== Cell A Collection ==================
# Create a collection for cell A
goo.make_collection("A_Cells", type = 'cell')
# Define cell A1
goo.make_cell("cell_A1", loc=(-1,0,0), collection="A_Cells")
# Define cell A2
goo.make_cell("cell_A2", loc=(1,0,0), collection="A_Cells")
# Define cell A3
goo.make_cell("cell_A3", loc=(0,-1,1.2), collection="A_Cells")
# Define cell A4
goo.make_cell("cell_A4", loc=(0,1,1.2), collection="A_Cells")

#================== Force A Collection ==================
# Create a collection for force A
goo.make_collection("A_Forces", type = 'force')
# Define force A1
fA1 = goo.make_force("force_A1", "cell_A1", strength=-1850, falloff=1, collection="A_Forces")
# Define force A2
fA2 = goo.make_force("force_A2", "cell_A2", strength=-1850, falloff=1, collection="A_Forces")
# Define force A3
fA3 = goo.make_force("force_A3", "cell_A3", strength=-1850, falloff=1, collection="A_Forces")
# Define force A4
fA4 = goo.make_force("force_A4", "cell_A4", strength=-1850, falloff=1, collection="A_Forces")

#================== Simulation setup ==================
handlers = goo.handler_class()
handlers.launch_simulation(start = 1, 
                           end = 500, 
                           filepath = "//your_file_path//", 
                           data = False, 
                           adhesion = True, 
                           growth = False, 
                           division = False, 
                           motility = False
                           )