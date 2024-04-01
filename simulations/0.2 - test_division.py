from importlib import reload
import goo
from goo.division import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

cell = goo.create_cell("cell", (1, 2, 3), subdivisions=3, physics_on=False)
print(f"Initial volume: {cell.volume()}")
logic = BisectDivisionLogic()
mother, daughter = cell.divide(logic)
logic.flush()
print(mother.name, daughter.name)
print(f"Volumes after division: {mother.volume()}, {daughter.volume()}")
