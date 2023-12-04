from goo import goo
from importlib import reload

reload(goo)
goo.setup_world(seed=1)

# Cells
# Define cell A1
goo.make_cell("cell_A1", loc=(0, 10, 0), type="cellsA")
# Define cell A2
goo.make_cell("cell_A2", loc=(10, 10, 0), type="cellsA")
# Define cell A1
goo.make_cell("cell_A3", loc=(0, 0, 0), type="cellsA")
# Define cell A2
goo.make_cell("cell_A4", loc=(10, 0, 0), type="cellsA")

# Forces

# Define forces A1
goo.add_homo_adhesion("cell_A1", -2000)
goo.add_motion('cell_A1', -0, distribution='uniform', size=0.1)
# Define forces A2
goo.add_homo_adhesion("cell_A2", -2000)
goo.add_motion('cell_A2', -500)
# Define forces A1
goo.add_homo_adhesion("cell_A3", -2000)
goo.add_motion('cell_A3', -1000)
# Define forces A2
goo.add_homo_adhesion("cell_A4", -2000)
goo.add_motion('cell_A4', -2000)

# Simulation setup
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=500,  # default, 250
                           filepath="C:\\tmp\\differential_motion\\data2", 
                           adhesion=True,  # default, True
                           data=True,  # default, False
                           growth=True, 
                           division=False, 
                           motility=True
                           )
