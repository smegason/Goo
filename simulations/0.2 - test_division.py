from importlib import reload
import goo
from goo.cell import *
from goo.reloader import *
from goo.division import *

reload(goo)
reset_modules()
reset_scene()

cell = create_cell("cell", (1, 2, 3), subdivisions=3)
print(f"Initial volume: {cell.volume()}")
mother, daughter = cell.divide(BisectDivisionLogic())
print(mother.name, daughter.name)
print(f"Volumes after division: {mother.volume()}, {daughter.volume()}")
