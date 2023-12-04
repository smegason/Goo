from goo import goo
from importlib import reload
reload(goo)
goo.setup_world(seed=1)

  
# Cells A and B
matB = ('red', 0.5, 0, 0)

goo.make_cell("cell_A1", loc=(-3, 1.5, 1.5), type="cellsA")
goo.make_cell("cell_A2", loc=(-3, -1.5, 1.5), type="cellsB", material=matB)
goo.make_cell("cell_A3", loc=(0, 1.5, 1.5), type="cellsB", material=matB)
goo.make_cell("cell_A4", loc=(0, -1.5, 1.5), type="cellsA")
goo.make_cell("cell_A5", loc=(3, 1.5, 1.5), type="cellsA")
goo.make_cell("cell_A6", loc=(3, -1.5, 1.5), type="cellsB", material=matB)
goo.make_cell("cell_A7",  loc=(-3, 1.5, -1.5), type="cellsB", material=matB)
goo.make_cell("cell_A8", loc=(-3, -1.5, -1.5), type="cellsA")
goo.make_cell("cell_A9", loc=(0, 1.5, -1.5), type="cellsA")
goo.make_cell("cell_A10", loc=(0, -1.5, -1.5), type="cellsB", material=matB)
goo.make_cell("cell_A11", loc=(3, 1.5, -1.5), type="cellsB", material=matB)
goo.make_cell("cell_A12", loc=(3, -1.5, -1.5), type="cellsA")
goo.make_cell("cell_A13", loc=(3, 4.5, 1.5), type="cellsB", material=matB)
goo.make_cell("cell_A14", loc=(3, 4.5, -1.5), type="cellsA")
goo.make_cell("cell_A15", loc=(0, 4.5, 1.5), type="cellsA")
goo.make_cell("cell_A16", loc=(0, 4.5, -1.5), type="cellsB", material=matB)
goo.make_cell("cell_A17", loc=(-3, 4.5, 1.5), type="cellsB", material=matB)
goo.make_cell("cell_A18",  loc=(-3, 4.5, -1.5), type="cellsA")

# Forces A and B
 
homoA = 2000
motion = 500
size = 1

goo.add_homo_adhesion('cell_A1', -homoA)
goo.add_motion('cell_A1', -motion, size=size)

goo.add_homo_adhesion('cell_A2', -homoA)
goo.add_motion('cell_A2', -motion, size=size)

goo.add_homo_adhesion('cell_A3', -homoA)
goo.add_motion('cell_A3', -motion, size=size)

goo.add_homo_adhesion('cell_A4', -homoA)
goo.add_motion('cell_A4', -motion, size=size)

goo.add_homo_adhesion('cell_A5', -homoA)
goo.add_motion('cell_A5', -motion, size=size)

goo.add_homo_adhesion('cell_A6', -homoA)
goo.add_motion('cell_A6', -motion, size=size)

goo.add_homo_adhesion('cell_A7', -homoA)
goo.add_motion('cell_A7', -motion, size=size)

goo.add_homo_adhesion('cell_A8', -homoA)
goo.add_motion('cell_A8', -motion, size=size)

goo.add_homo_adhesion('cell_A9', -homoA)
goo.add_motion('cell_A9', -motion, size=size)

goo.add_homo_adhesion('cell_A10', -homoA)
goo.add_motion('cell_A10', -motion, size=size)

goo.add_homo_adhesion('cell_A11', -homoA)
goo.add_motion('cell_A11', -motion, size=size)

goo.add_homo_adhesion('cell_A12', -homoA)
goo.add_motion('cell_A12', -motion, size=size)

goo.add_homo_adhesion('cell_A13', -homoA)
goo.add_motion('cell_A13', -motion, size=size)

goo.add_homo_adhesion('cell_A14', -homoA)
goo.add_motion('cell_A14', -motion, size=size)

goo.add_homo_adhesion('cell_A15', -homoA)
goo.add_motion('cell_A15', -motion, size=size)

goo.add_homo_adhesion('cell_A16', -homoA)
goo.add_motion('cell_A16', -motion, size=size)

goo.add_homo_adhesion('cell_A17', -homoA)
goo.add_motion('cell_A17', -motion, size=size)

goo.add_homo_adhesion('cell_A18', -homoA)
goo.add_motion('cell_A18', -motion, size=size)

# Reflective boundaries

goo.add_boundaries(loc=(0, 1.5, 0), size=(6, 6, 4))

# Simulation setup
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=1000,  # default, 250
                           filepath="/Users/antoine/Harvard/MegasonLab/presentations/20231205_ASCB_RikkiGarner/20231129_cellsorting_distr0.1/data1", 
                           adhesion=True,  # default, True
                           data=True,  # default, False
                           growth=True, 
                           motility=True
                           )
