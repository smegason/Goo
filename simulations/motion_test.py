from goo import goo
from importlib import reload

reload(goo)
goo.setup_world(seed=1)

# Define cell A1
goo.make_cell("cell_A1", loc=(0, 0, 0), type="cellsA")

# Define forces A1
goo.add_homo_adhesion("cell_A1", -5000)
goo.add_motion('cell_A1', -2500, distribution='uniform')

# Simulation setup
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=200,  # default, 250
                           filepath="/Users/antoine/Harvard/MegasonLab/GPU_backup/AntoineRuzette/goo/data/motion_diffusion/20240402_motion_1/20240402_motion_1", 
                           adhesion=True,  # default, True
                           data=True,  # default, False
                           growth='tissue', 
                           division=False, 
                           motility=True,
                           target_volume=20
                           )
