# hello_world.py - creates one cell

from goo import goo
from importlib import reload
reload(goo)

goo.setup_world()

goo.add_material_cell("CellGreen", 0.007, 0.300, 0.005)
goo.add_material_cell("CellBlue", 0.007, 0.021, 0.300)

cell = goo.Cell(name_string="Cell1_", loc=(0, 0, 0), material="CellGreen")
goo.make_cell(cell)

cell2 = goo.Cell(name_string="Cell2_", loc=(2, 0, 0), material="CellBlue")
goo.make_cell(cell2)
