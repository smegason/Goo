from importlib import reload
import goo
from goo.division import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

cell = goo.create_cell(
    "cell",
    (1, 2, 3),
    subdivisions=3,
    scale=(2, 3, 1),
    physics_on=False,
)
print(f"Initial volume: {cell.get_volume()}")
logic = BisectDivisionLogic(margin=0.025)
# logic = BooleanDivisionLogic()
mother, daughter = cell.divide(logic)
logic.flush()

print(mother.name, daughter.name)
print(f"Volumes after division: {mother.get_volume()}, {daughter.get_volume()}")
