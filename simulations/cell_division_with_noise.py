from goo import goo
from importlib import reload
reload(goo)
goo.setup_world()

# Define cell A1
goo.make_cell("cell_A1", loc=(0, 0, 0), type="cellsA", scale=(1.67, 1.67, 1.67))

# Simulation setup 
handlers = goo.handler_class()
handlers.launch_simulation(start=1,
                           end=130, 
                           filepath="/Users/antoine/Harvard/MegasonLab/GPU_backup/AntoineRuzette/goo/data/division/20240304_desynchronized_division/20240304_desynchronized_division_2", 
                           adhesion=False,
                           data=True,
                           growth=True, 
                           motility=False, 
                           division=True, 
                           target_volume=30, 
                           division_type='time', 
                           cell_cycle_time=20, 
                           growth_rate=1, 
                           dt_physics=1, 
                           cell_cycle_var=0.001
                           )
