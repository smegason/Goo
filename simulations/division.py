from goo import goo
from importlib import reload
reload(goo)
goo.setup_world(seed=1)

# Cell A
goo.make_cell("cell_A1", loc=(0, 0, 0), type="cellsA")
 
# Force A
homoA = 2000
motion = 500
size = 1
goo.add_homo_adhesion('cell_A1', -homoA)


# Simulation setup
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=10000,  # default, 250
                           filepath="/Users/antoine/Harvard/MegasonLab/presentations/20231206_soma/division_example.mp4", 
                           adhesion=True,  # default, True
                           data=True,  # default, False
                           growth=True, 
                           motility=False, 
                           division=True
                           )

