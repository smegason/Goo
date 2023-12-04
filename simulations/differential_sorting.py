from goo import goo
from importlib import reload

reload(goo)
goo.setup_world(seed=1)

# Cells A 

# Define cell A1
goo.make_cell("cell_A1", loc=(2, 0, 0), type="cellsA")
# Define cell A2
goo.make_cell("cell_A2", loc=(0, 2, 0), type="cellsA")
# Define cell A3
goo.make_cell("cell_A3", loc=(4, 2, 0), type="cellsA")
# Define cell A4
goo.make_cell("cell_A4", loc=(0, -2, 0), type="cellsA")
# Define cell A5
goo.make_cell("cell_A5", loc=(4, -2, 0), type="cellsA")

# Define cell A6
goo.make_cell("cell_A6", loc=(0, 0, 2), type="cellsA")
# Define cell A7
goo.make_cell("cell_A7", loc=(4, 0, 2), type="cellsA")
# Define cell A8
goo.make_cell("cell_A8", loc=(2, 2, 2), type="cellsA")
# Define cell A9
goo.make_cell("cell_A9", loc=(2, -2, 2), type="cellsA")

# Define cell A10
goo.make_cell("cell_A10", loc=(0, 0, -2), type="cellsA")
# Define cell A11
goo.make_cell("cell_A11", loc=(4, 0, -2), type="cellsA")
# Define cell A12
goo.make_cell("cell_A12", loc=(2, 2, -2), type="cellsA")
# Define cell A13
goo.make_cell("cell_A13", loc=(2, -2, -2), type="cellsA")

# Cells B

# Define cell B1
goo.make_cell("cell_B1", loc=(0, 0, 0), material=('red', 0.2, 0, 0), type="cellsB")
# Define cell B2
goo.make_cell("cell_B2", loc=(4, 0, 0), material=('red', 0.2, 0, 0), type="cellsB")
# Define cell B3
goo.make_cell("cell_B3", loc=(2, 2, 0), material=('red', 0.2, 0, 0), type="cellsB")
# Define cell B4
goo.make_cell("cell_B4", loc=(2, -2, 0), material=('red', 0.2, 0, 0), type="cellsB")
# Define cell B5
goo.make_cell("cell_B5", loc=(2, 0, -2), material=('red', 0.2, 0, 0), type="cellsB")
# Define cell B6
goo.make_cell("cell_B6", loc=(0, 2, -2), material=('red', 0.2, 0, 0), type="cellsB")
# Define cell B7
goo.make_cell("cell_B7", loc=(4, 2, -2), material=('red', 0.2, 0, 0), type="cellsB")
# Define cell B8
goo.make_cell("cell_B8", loc=(0, -2, -2), material=('red', 0.2, 0, 0), type="cellsB")
# Define cell B9
goo.make_cell("cell_B9", loc=(4, -2, -2), material=('red', 0.2, 0, 0), type="cellsB")
# Define cell B10
goo.make_cell("cell_B10", loc=(2, 0, 2), material=('red', 0.2, 0, 0), type="cellsB")
# Define cell B11
goo.make_cell("cell_B11", loc=(0, 2, 2), material=('red', 0.2, 0, 0), type="cellsB")
# Define cell B12
goo.make_cell("cell_B12", loc=(4, 2, 2), material=('red', 0.2, 0, 0), type="cellsB")
# Define cell B13
goo.make_cell("cell_B13", loc=(0, -2, 2), material=('red', 0.2, 0, 0), type="cellsB")
# Define cell B14
goo.make_cell("cell_B14", loc=(4, -2, 2), material=('red', 0.2, 0, 0), type="cellsB")

# Forces A 

homoA = 2500
homoB = 2500
heteroAB = 500
motion = 500

# Define force A1
goo.add_homo_adhesion('cell_A1', -homoA)
goo.add_hetero_adhesion('cell_A1', 'cellsB', -heteroAB)
goo.add_motion('cell_A1', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A2', -homoA)
goo.add_hetero_adhesion('cell_A2', 'cellsB', -heteroAB)
goo.add_motion('cell_A2', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A3', -homoA)
goo.add_hetero_adhesion('cell_A3', 'cellsB', -heteroAB)
goo.add_motion('cell_A3', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A4', -homoA)
goo.add_hetero_adhesion('cell_A4', 'cellsB', -heteroAB)
goo.add_motion('cell_A4', -motion)
# Define force A1
goo.add_homo_adhesion('cell_A5', -homoB)
goo.add_hetero_adhesion('cell_A5', 'cellsB', -heteroAB)
goo.add_motion('cell_A5', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A6', -homoB)
goo.add_hetero_adhesion('cell_A6', 'cellsB', -heteroAB)
goo.add_motion('cell_A6', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A7', -homoB)
goo.add_hetero_adhesion('cell_A7', 'cellsB', -heteroAB)
goo.add_motion('cell_A7', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A8', -homoB)
goo.add_hetero_adhesion('cell_A8', 'cellsB', -heteroAB)
goo.add_motion('cell_A8', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A9', -homoA)
goo.add_hetero_adhesion('cell_A9', 'cellsB', -heteroAB)
goo.add_motion('cell_A9', -motion)
# Define force A1
goo.add_homo_adhesion('cell_A10', -homoB)
goo.add_hetero_adhesion('cell_A10', 'cellsB', -heteroAB)
goo.add_motion('cell_A10', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A11', -homoB)
goo.add_hetero_adhesion('cell_A11', 'cellsB', -heteroAB)
goo.add_motion('cell_A11', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A12', -homoB)
goo.add_hetero_adhesion('cell_A12', 'cellsB', -heteroAB)
goo.add_motion('cell_A12', -motion)
# Define force A2
goo.add_homo_adhesion('cell_A13', -homoB)
goo.add_hetero_adhesion('cell_A13', 'cellsB', -heteroAB)
goo.add_motion('cell_A13', -motion)

# Forces B 

# Define force B1
goo.add_homo_adhesion('cell_B1', -homoB)
goo.add_hetero_adhesion('cell_B1', 'cellsA', -heteroAB)
goo.add_motion('cell_B1', -motion)
# Define force B2
goo.add_homo_adhesion('cell_B2', -homoB)
goo.add_hetero_adhesion('cell_B2', 'cellsA', -heteroAB)
goo.add_motion('cell_B2', -motion)
# Define force B3
goo.add_homo_adhesion('cell_B3', -homoB)
goo.add_hetero_adhesion('cell_B3', 'cellsA', -heteroAB)
goo.add_motion('cell_B3', -motion)
# Define force B4
goo.add_homo_adhesion('cell_B4', -homoB)
goo.add_hetero_adhesion('cell_B4', 'cellsA', -heteroAB)
goo.add_motion('cell_B4', -motion)
# Define force B5
goo.add_homo_adhesion('cell_B5', -homoB)
goo.add_hetero_adhesion('cell_B5', 'cellsA', -heteroAB)
goo.add_motion('cell_B5', -motion)
# Define force B6
goo.add_homo_adhesion('cell_B6', -homoB)
goo.add_hetero_adhesion('cell_B6', 'cellsA', -heteroAB)
goo.add_motion('cell_B6', -motion)
# Define force B7
goo.add_homo_adhesion('cell_B7', -homoB)
goo.add_hetero_adhesion('cell_B7', 'cellsA', -heteroAB)
goo.add_motion('cell_B7', -motion)
# Define force B8
goo.add_homo_adhesion('cell_B8', -homoB)
goo.add_hetero_adhesion('cell_B8', 'cellsA', -heteroAB)
goo.add_motion('cell_B8', -motion)
# Define force B9
goo.add_homo_adhesion('cell_B9', -homoB)
goo.add_hetero_adhesion('cell_B9', 'cellsA', -heteroAB)
goo.add_motion('cell_B9', -motion)
# Define force B10
goo.add_homo_adhesion('cell_B10', -homoB)
goo.add_hetero_adhesion('cell_B10', 'cellsA', -heteroAB)
goo.add_motion('cell_B10', -motion)
# Define force B11
goo.add_homo_adhesion('cell_B11', -homoB)
goo.add_hetero_adhesion('cell_B11', 'cellsA', -heteroAB)
goo.add_motion('cell_B11', -motion)
# Define force B12
goo.add_homo_adhesion('cell_B12', -homoB)
goo.add_hetero_adhesion('cell_B12', 'cellsA', -heteroAB)
goo.add_motion('cell_B12', -motion)
# Define force B13
goo.add_homo_adhesion('cell_B13', -homoB)
goo.add_hetero_adhesion('cell_B13', 'cellsA', -heteroAB)
goo.add_motion('cell_B13', -motion)
# Define force B14
goo.add_homo_adhesion('cell_B14', -homoB)
goo.add_hetero_adhesion('cell_B14', 'cellsA', -heteroAB)
goo.add_motion('cell_B14', -motion)

# Simulation setup 
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=500,  # default, 250
                           filepath="C:\\Users\\anr9744\\Projects\\Goo\\data\\cell_sorting", 
                           adhesion=True,  # default, True
                           data=False,  # default, False
                           growth=True,
                           motility=True,
                           )
