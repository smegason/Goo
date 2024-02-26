from goo import goo
from importlib import reload
# import bpy
reload(goo)
goo.setup_world(seed=1)

# Cell A
goo.make_cell("cell_A1", loc=(-5, 2.5, 0), type="cellsA")
goo.make_cell("cell_A2", loc=(0, 2.5, 0), type="cellsA")
goo.make_cell("cell_A3", loc=(5, 2.5, 0), type="cellsA")
goo.make_cell("cell_A4", loc=(-5, -2.5, 0), type="cellsA")
goo.make_cell("cell_A5", loc=(0, -2.5, 0), type="cellsA")
goo.make_cell("cell_A6", loc=(5, -2.5, 0), type="cellsA")

goo.add_turbulence_motion(strength=-1000, seed=100)

goo.add_sphere_boundaries(loc=(0, 0, 0), radius=10)

# Simulation setup
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=500,  # default, 250
                           filepath="/Users/antoine/Harvard/MegasonLab/GPU_backup/AntoineRuzette/goo/data/motion_diffusion/20240211_turbulencemotioninsphere_seeds/20240211_turbulencemotioninsphere_seed100", 
                           adhesion=False,  # default, True
                           data=True,  # default, False
                           growth=True, 
                           motility=False, 
                           division=False, 
                           target_volume=50,
                           growth_rate=1, 
                           growth_type='linear'
                           )
