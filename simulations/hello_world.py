# hello_world.py - creates one cell
from goo import goo

goo.setup_world()

cell = goo.Cell(name_string="Cell1_", loc=(0, 0, 0), material="CellGreen")
goo.make_cell(cell)
