# Goal: test the ability to import and reload modules

from importlib import reload
import goo

reload(goo)
goo.reset_modules()
goo.reset_scene()

celltype = goo.CellType("A", pattern="simple")
cell = celltype.create_cell("cell", (0, 0, 0), size=1)

print(cell.area())
print(cell.volume())
