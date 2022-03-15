# hello_world.py - creates one cell
from Goo import Goo

Goo.setup_world()

cell = Goo.Cell(name_string = "Cell1_", loc = (0, 0, 0), material = "CellGreen")
Goo.make_cell(cell)
