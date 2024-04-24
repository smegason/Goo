from goo import goo
from importlib import reload
reload(goo)
goo.setup_world()

# Cells A 
goo.make_cell("cell_A1", loc=(-1.9, 0, 0), type="A_Cells", scale=(1.67, 1.67, 1.67))
goo.make_cell("cell_A2", loc=(1.9, 0, 0), type="A_Cells", scale=(1.67, 1.67, 1.67))

# Forces A
homoA = 7500

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)
goo.add_homo_adhesion('cell_A2', -homoA)

# Simulation setup 
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=110,  # default, 250
                           filepath="", 
                           adhesion=True,  # default, True
                           data=True,  # default, False
                           growth=True, 
                           motility=False, 
                           division=True, 
                           target_volume=30, 
                           cell_cycle_time=30,
                           cell_cycle_var=0
                           )
