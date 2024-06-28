from importlib import reload
import goo
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()


celltype = goo.CellType("A")
cell = celltype.create_cell("cellA", (0, 0, 0))
cell.recolor((1, 0, 0))
