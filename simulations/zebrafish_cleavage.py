# zebrafish_cleavage.py - simulates the first 12 cell cycles of zebrafish development
import bpy
from goo import goo
from importlib import reload
reload(goo)

goo.setup_world()

# make the blastomere
cell = goo.Cell("cell_", loc=(0, 0, 0))

# make the yolk
yolk = goo.Cell("yolk", loc=(0, 0, -1))

# make cell and yolk
goo.make_cell(cell)
goo.make_cell(yolk)
