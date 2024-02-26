from goo import goo
from importlib import reload
# import bpy
reload(goo)
goo.setup_world(seed=1)

# Cell A
goo.make_cell("cell_A1", loc=(-5, -3, -3), type="cellsA")
goo.make_cell("cell_A2", loc=(0, -3, -3), type="cellsA")
goo.make_cell("cell_A3", loc=(5, -3, -3), type="cellsA")
goo.make_cell("cell_A4", loc=(-5, -3, 3), type="cellsA")
goo.make_cell("cell_A5", loc=(0, -3, 3), type="cellsA")
goo.make_cell("cell_A6", loc=(5, -3, 3), type="cellsA")
goo.make_cell("cell_A7", loc=(-5, 3, -3), type="cellsA")
goo.make_cell("cell_A8", loc=(0, 3, -3), type="cellsA")
goo.make_cell("cell_A9", loc=(5, 3, -3), type="cellsA")
goo.make_cell("cell_A10", loc=(-5, 3, 3), type="cellsA")
goo.make_cell("cell_A11", loc=(0, 3, 3), type="cellsA")
goo.make_cell("cell_A12", loc=(5, 3, 3), type="cellsA")

# Force A
homoA = 2000
goo.add_homo_adhesion('cell_A1', -homoA)
goo.add_homo_adhesion('cell_A2', -homoA)
goo.add_homo_adhesion('cell_A3', -homoA)
goo.add_homo_adhesion('cell_A4', -homoA)
goo.add_homo_adhesion('cell_A5', -homoA)
goo.add_homo_adhesion('cell_A6', -homoA)
goo.add_homo_adhesion('cell_A7', -homoA)
goo.add_homo_adhesion('cell_A8', -homoA)
goo.add_homo_adhesion('cell_A9', -homoA)
goo.add_homo_adhesion('cell_A10', -homoA)
goo.add_homo_adhesion('cell_A11', -homoA)
goo.add_homo_adhesion('cell_A12', -homoA)

goo.add_turbulence_motion(strength=-20000)

goo.add_sphere_boundaries(loc=(0, 0, 0), radius=10)

# Simulation setup
handlers = goo.handler_class()
handlers.launch_simulation(start=1,  # default, 1
                           end=10000,  # default, 250
                           filepath="/Users/antoine/Harvard/MegasonLab/GPU_backup/AntoineRuzette/goo/data/division/20240102_division_tests/test01", 
                           adhesion=True,  # default, True
                           data=True,  # default, False
                           growth=True, 
                           motility=False, 
                           division=False, 
                           target_volume=50,
                           growth_rate=1, 
                           growth_type='linear'
                           )
