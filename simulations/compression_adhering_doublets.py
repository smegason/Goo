from goo import goo
from importlib import reload

reload(goo)
goo.setup_world()

# Cells A 
goo.make_cell("cell_A1", loc=(-1, 0, 0), type="cellsA")
goo.make_cell("cell_A2", loc=(1, 0, 0), type="cellsA")

goo.make_cell("cell_B3", loc=(-1, 2.2, 0), material=('green', 0, .1, 0), type="cellsB")
goo.make_cell("cell_B4", loc=(1, 2.2, 0), material=('green', 0, .1, 0), type="cellsB")

goo.make_cell("cell_C5", loc=(-1, 4.4, 0), material=('red', .1, 0, 0), type="cellsC")
goo.make_cell("cell_C6", loc=(1, 4.4, 0), material=('red', .1, 0, 0), type="cellsC")

# Forces 

homoA = 1000
homoB = 2000
homoC = 3000

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)
goo.add_homo_adhesion('cell_A2', -homoA)

goo.add_homo_adhesion('cell_B3', -homoB)
goo.add_homo_adhesion('cell_B4', -homoB)

goo.add_homo_adhesion('cell_C5', -homoC)
goo.add_homo_adhesion('cell_C6', -homoC)

# Simulation setup
handlers = goo.handler_class()
handlers.launch_simulation(start=1,
                           end=500,
                           filepath="/Users/antoine/Harvard/MegasonLab/GPU_backup/AntoineRuzette/goo/data/3D-physics/compression/20240122_compression_adhering_doublets/20240122_compression_adhering_doublets", 
                           adhesion=True,
                           data=True,
                           growth=True, 
                           division=False, 
                           motility=False, 
                           volume_scale=1
                           )
