from goo import goo
from importlib import reload

reload(goo)
goo.setup_world()

# Cells A 
goo.make_cell("cell_A1", loc=(-1.65, 0, 0), type="cellsA", scale=(1.67, 1.67, 1.67))
goo.make_cell("cell_A2", loc=(1.65, 0, 0), type="cellsA", scale=(1.67, 1.67, 1.67))
# Forces 

homoA = 10000

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)
goo.add_homo_adhesion('cell_A2', -homoA)


# Simulation setup
handlers = goo.handler_class()
handlers.launch_simulation(start=1,
                           end=500,
                           filepath="/Users/antoine/Harvard/MegasonLab/GPU_backup/AntoineRuzette/goo/data/3D-physics/compression/20240122_compression_adhering_doublets/20240122_compression_adhering_doublets", 
                           adhesion=True,
                           data=True,
                           growth='tissue', 
                           division=False, 
                           motility=False, 
                           target_volume=30
                           )
