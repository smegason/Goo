from goo import goo
from importlib import reload
reload(goo)
goo.setup_world()

  
# Cells A 
# Define cell A1
goo.make_cell("cell_A1", 
              loc=(0, 0, 0), 
              type="cellsA", 
              scale=(2, 2, 2),
              stiffness=1)

goo.make_yolk("yolk", 
              loc=(0, 0, -5.75), 
              type="cellsA", 
              scale=(0.38, 0.38, 0.38), 
              adhesion_strength=12000,
              stiffness=1)

# Forces A
homoA = 10000

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)

# Simulation setup 
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=70,  # default, 250
                           filepath="/Users/antoine/Harvard/MegasonLab/GPU_backup/AntoineRuzette/goo/data/division/20240318_cleaving_embryo_synchronous/20240318_cleaving_embryo_synchronous_adh10k", 
                           adhesion=True,  # default, True
                           data=True,  # default, False
                           growth=True, 
                           motility=False, 
                           division=True, 
                           target_volume=50, 
                           cell_cycle_time=30,
                           cell_cycle_var=0
                           )
