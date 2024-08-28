# Goal: test the ability to perform divisions,
# and that the division is done correctly.

from importlib import reload
import goo
from goo.division import *
from math import pi

reload(goo)
goo.reset_modules()
goo.reset_scene()

rot = (pi / 4, pi / 4, pi / 4)
cell = goo.create_cell(
    "cell",
    (0, 0, 0),
    subdivisions=4,
    scale=(2, 1, 3),
    rotation=rot,
    physics_on=False,
)

print(f"Initial volume: {cell.volume()}")
logic = BisectDivisionLogic(margin=0.025)
# logic = BooleanDivisionLogic()
mother, daughter = cell.divide(logic)
logic.flush()

print(mother.name, daughter.name)
print(f"Volumes after division: {mother.volume()}, {daughter.volume()}")
