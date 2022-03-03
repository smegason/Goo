# hello_world.py - creates one cell

from Goo import Goo

Goo.setup_world()
cell = Goo.Cell(name_string = "Cell_", loc = (0, 0, 0))
Goo.make_cell(cell)