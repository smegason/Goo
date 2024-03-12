from goo import goo
from importlib import reload
reload(goo)
goo.setup_world()

  
# Cells A 
# Define cell A1
goo.make_cell("cell_A1", loc=(0, 0, 0), type="A_Cells", scale=(1.67, 1.67, 1.67))

# Simulation setup 
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=500,  # default, 250
                           filepath="C:\\Users\\tmp\\", 
                           adhesion=False,  # default, True
                           data=True,  # default, False
                           growth=True, 
                           motility=False, 
                           division=True, 
                           target_volume=30
                           )
