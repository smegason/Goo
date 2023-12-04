# doublet.py - simulates 2 cells adhering to each other
# with different balances of cortical tension and cell adhesion

from goo import goo
from importlib import reload

reload(goo)
goo.setup_world()


# Cells A

# Define cell A1
goo.make_cell("cell_A1", loc=(-1.36, -1, 0.2), type="cellsA")
# Define cell A2
goo.make_cell("cell_A2", loc=(-0.2, -1, -1.36), type="cellsA")
# Define cell A3
goo.make_cell("cell_A3", loc=(1.36, -1, -0.2), type="cellsA")
# Define cell A4
goo.make_cell("cell_A4", loc=(0.2, -1, 1.36), type="cellsA")
# Define cell A5
goo.make_cell("cell_A5", loc=(-0.73, 1, 1.1), type="cellsA")
# Define cell A6
goo.make_cell("cell_A6", loc=(1.2, 1, 0.84), type="cellsA")
# Define cell A7
goo.make_cell("cell_A7", loc=(0.92, 1, -1.08), type="cellsA")
# Define cell A8
goo.make_cell("cell_A8", loc=(-1, 1, -0.82), type="cellsA")


# Forces A

homoA = 2000
motion = 500

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)
goo.add_motion('cell_A1', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A2', -homoA)
goo.add_motion('cell_A2', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A3', -homoA)
goo.add_motion('cell_A3', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A4', -homoA)
goo.add_motion('cell_A4', -motion)


# Define force A1
goo.add_homo_adhesion('cell_A5', -homoA)
goo.add_motion('cell_A5', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A6', -homoA)
goo.add_motion('cell_A6', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A7', -homoA)
goo.add_motion('cell_A7', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A8', -homoA)
goo.add_motion('cell_A8', -motion)


# Simulation setup
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=500,  # default, 250
                           filepath="C:\\tmp\\sorting_test_rendering\\data5", 
                           adhesion=True,  # default, True
                           data=False,  # default, False
                           growth=True, 
                           division=False, 
                           motility=True
                           )
