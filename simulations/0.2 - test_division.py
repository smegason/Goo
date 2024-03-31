from importlib import reload
import goo
from goo.cell import *
from goo.reloader import *
from goo.division import *
from goo.force import *

reload(goo)
reset_modules()
reset_scene()

cell = create_cell("cell", (1, 2, 3), subdivisions=3, physics_on=False)
print(f"Initial volume: {cell.volume()}")
logic = BisectDivisionLogic()
mother, daughter = cell.divide(logic)
logic.flush()
print(mother.name, daughter.name)
print(f"Volumes after division: {mother.volume()}, {daughter.volume()}")
