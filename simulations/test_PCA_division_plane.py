from goo import goo
from importlib import reload
reload(goo)
goo.setup_world()

  
# Cells A 
# Define cell A1
goo.make_cell("cell_A1", loc=(0, 0, 0), type="A_Cells", 
              scale=(1, 0.8, 0.8), rotation=(0, 45, 0))
goo.make_cell("cell_A2", loc=(0, 4, 0), type="A_Cells", 
              scale=(0.8, 0.8, 1), rotation=(45, 0, 0))
goo.make_cell("cell_A3", loc=(0, -4, 0), type="A_Cells", 
              scale=(0.8, 1, 0.8), rotation=(0, 0, 45))

# Simulation setup 
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=500,  # default, 250
                           filepath="C:\\Users\\tmp\\", 
                           adhesion=True,  # default, True
                           data=False,  # default, False
                           growth=True, 
                           motility=False, 
                           division=False, 
                           test=True
                           )
