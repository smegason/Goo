from goo import goo
from importlib import reload

reload(goo)
goo.setup_world()

# Cells A 
# Define cell A1
goo.make_cell("cell_A1", loc=(-1.75, 0, 0), type="cellsA", scale=(1.6, 1.6, 1.6))
# Define cell A2
goo.make_cell("cell_A2", loc=(1.75, 0, 0), type="cellsA", scale=(1.6, 1.6, 1.6))

# Define cell B1
goo.make_cell("cell_B3", loc=(-1.75, 2.2, 0), material=('green', 0, .1, 0), type="cellsB", scale=(1.6, 1.6, 1.6))
# Define cell B2
goo.make_cell("cell_B4", loc=(1.75, 2.2, 0), material=('green', 0, .1, 0), type="cellsB", scale=(1.6, 1.6, 1.6))

# Define cell C1
goo.make_cell("cell_C5", loc=(-1.75, 4.4, 0), material=('red', .1, 0, 0), type="cellsC", scale=(1.6, 1.6, 1.6))
# Define cell C2
goo.make_cell("cell_C6", loc=(1.75, 4.4, 0), material=('red', .1, 0, 0), type="cellsC", scale=(1.6, 1.6, 1.6))

# Forces 

homoA = 5000
homoB = 7500
homoC = 10000

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)
# goo.add_hetero_adhesion('cell_A1', 'cellsB', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_A2', -homoA)
# goo.add_hetero_adhesion('cell_A2', 'cellsB', -heteroAB)
# Define force B1
goo.add_homo_adhesion('cell_B3', -homoB)
# goo.add_hetero_adhesion('cell_A3', 'cellsB', -heteroAB)
# Define force B2
goo.add_homo_adhesion('cell_B4', -homoB)
# goo.add_hetero_adhesion('cell_A4', 'cellsB', -heteroAB)
# Define force C1
goo.add_homo_adhesion('cell_C5', -homoC)
# goo.add_hetero_adhesion('cell_B5', 'cellsA', -heteroAB)
# Define force C2
goo.add_homo_adhesion('cell_C6', -homoC)
# goo.add_hetero_adhesion('cell_B6', 'cellsA', -heteroAB)


# Simulation setup
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=500,  # default, 250
                           filepath="C:\\tmp\\sorting_test_rendering\\data5", 
                           adhesion=True,  # default, True
                           data=False,  # default, False
                           growth=True, 
                           division=False, 
                           motility=False, 
                           target_volume=20
                           )
