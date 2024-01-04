from goo import goo
from importlib import reload

reload(goo)
goo.setup_world(seed=1)

# Cells A  
# Define cell A1
goo.make_cell("cell_A1", loc=(-1, 0, 0), type="cellsA")
# Define cell A2
goo.make_cell("cell_A2", loc=(-3, 0, 0), type="cellsA")
# Define cell A3
goo.make_cell("cell_A3", loc=(-1, 2, 0), type="cellsA")
# Define cell A4
goo.make_cell("cell_A4", loc=(-3, 2, 0), type="cellsA")
# Define cell A1
goo.make_cell("cell_A5", loc=(-1, 0, -2), type="cellsA")
# Define cell A2
goo.make_cell("cell_A6", loc=(-3, 0, -2), type="cellsA")
# Define cell A3
goo.make_cell("cell_A7", loc=(-1, 2, -2), type="cellsA")
# Define cell A4
goo.make_cell("cell_A8", loc=(-3, 2, -2), type="cellsA")

# Define cell B1
goo.make_cell("cell_B1", loc=(1, 0, 0), material=('red', 0.1, 0, 0), type="cellsB")
# Define cell B2
goo.make_cell("cell_B2", loc=(3, 0, 0), material=('red', 0.1, 0, 0), type="cellsB")
# Define cell B3
goo.make_cell("cell_B3", loc=(1, 2, 0), material=('red', 0.1, 0, 0), type="cellsB")
# Define cell B4
goo.make_cell("cell_B4", loc=(3, 2, 0), material=('red', 0.1, 0, 0), type="cellsB")
# Define cell B1
goo.make_cell("cell_B5", loc=(1, 0, -2), material=('red', 0.1, 0, 0), type="cellsB")
# Define cell B2
goo.make_cell("cell_B6", loc=(3, 0, -2), material=('red', 0.1, 0, 0), type="cellsB")
# Define cell B3
goo.make_cell("cell_B7", loc=(1, 2, -2), material=('red', 0.1, 0, 0), type="cellsB")
# Define cell B4
goo.make_cell("cell_B8", loc=(3, 2, -2), material=('red', 0.1, 0, 0), type="cellsB")


# Forces A  

homoA = 2000
homoB = 2000
heteroAB = 4000

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)
goo.add_hetero_adhesion('cell_A1', 'cellsB', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_A2', -homoA)
goo.add_hetero_adhesion('cell_A2', 'cellsB', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_A3', -homoA)
goo.add_hetero_adhesion('cell_A3', 'cellsB', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_A4', -homoA)
goo.add_hetero_adhesion('cell_A4', 'cellsB', -heteroAB)

# Define force A1
goo.add_homo_adhesion('cell_B5', -homoB)
goo.add_hetero_adhesion('cell_B5', 'cellsA', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_B6', -homoB)
goo.add_hetero_adhesion('cell_B6', 'cellsA', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_B7', -homoB)
goo.add_hetero_adhesion('cell_B7', 'cellsA', -heteroAB)
# Define force A2
goo.add_homo_adhesion('cell_B8', -homoB)
goo.add_hetero_adhesion('cell_B8', 'cellsA', -heteroAB)


# Simulation setup 
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=300,  # default, 250
                           filepath="C:\\tmp\\adhesion_paper\\hetero_homo\\hetero", 
                           adhesion=True,  # default, True
                           data=True,  # default, False
                           growth=True, 
                           division=False, 
                           motility=False
                           )
