from goo import goo
from importlib import reload

reload(goo)
goo.setup_world()

# Cells A 
# Define cell A1
goo.make_cell("cell_A1", loc=(0, -1.85, 0), material=('blue', 0, 0, .5), type="cellsA", scale=(2, 2, 2))
# Define cell A2
goo.make_cell("cell_A2", loc=(0, 1.85, 0), material=('blue', 0, 0, .5), type="cellsA", scale=(2, 2, 2))

# Define cell B1
goo.make_cell("cell_B3", loc=(5, -1.85, 0), material=('green', 0, .5, 0), type="cellsB", scale=(2, 2, 2))
# Define cell B2
goo.make_cell("cell_B4", loc=(5, 1.85, 0), material=('green', 0, .5, 0), type="cellsB", scale=(2, 2, 2))

# Define cell C1
goo.make_cell("cell_C5", loc=(10, -1.85, 0), material=('red', .5, 0, 0), type="cellsC", scale=(2, 2, 2))
# Define cell C2
goo.make_cell("cell_C6", loc=(10, 1.85, 0), material=('red', .5, 0, 0), type="cellsC", scale=(2, 2, 2))

# Forces 

homoA = 2500
homoB = 5000
homoC = 7500

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)
# Define force A2
goo.add_homo_adhesion('cell_A2', -homoA)
# Define force B1
goo.add_homo_adhesion('cell_B3', -homoB)
# Define force B2
goo.add_homo_adhesion('cell_B4', -homoB)
# Define force C1
goo.add_homo_adhesion('cell_C5', -homoC)
# Define force C2
goo.add_homo_adhesion('cell_C6', -homoC)


# Simulation setup
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=200,  # default, 250
                           filepath="/Users/antoine/Harvard/MegasonLab/GPU_backup/AntoineRuzette/goo/data/adhesion/20240401_adhesion_doublets_4/20240401_adhesion_doublets", 
                           adhesion=True,  # default, True
                           data=False,  # default, False
                           growth='tissue', 
                           division=False, 
                           motility=False, 
                           target_volume=30, 
                           colorize='type'
                           )
