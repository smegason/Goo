# hello_world.py - creates one cell

import bpy
from Goo import Goo
from importlib import reload
reload(Goo)

Goo.setup_world()

Goo.add_material_cell("CellGreen",0.007,0.300,0.005)
Goo.add_material_cell("CellBlue",0.007, 0.021, 0.300)

cell = Goo.Cell(name_string = "Cell1_", loc = (0, 0, 0), material = "CellGreen")
Goo.make_cell(cell)

cell2 = Goo.Cell(name_string = "Cell2_", loc = (2, 0, 0), material = "CellBlue")
Goo.make_cell(cell2)