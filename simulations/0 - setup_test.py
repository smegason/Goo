# Goal: test the creation of a world, and population of cells within that world
# through base Cell creation and CellType creation.

from importlib import reload
import goo
from goo import Cell, CellType, reset

reload(goo)
reset()

cell = Cell("cell", (0, 0, 0))
