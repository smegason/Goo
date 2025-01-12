from importlib import reload
import goo
from goo import *

reload(goo)

reset_modules()
reset_scene()

sim = Simulator()

cellsA = CellType("A")
cellsA.homo_adhesion_strength = 2000

cellsB = CellType("B")
cellsB.homo_adhesion_strength = 2000

cellsA.set_hetero_adhesion(cellsB, 4000)

locsA = [
    (-1, 0, 0),
    (-3, 0, 0),
    (-1, 2, 0),
    (-3, 2, 0),
    (-1, 0, -2),
    (-3, 0, -2),
    (-1, 2, -2),
    (-3, 2, -2),
]

for i, loc in enumerate(locsA):
    cellsA.create_cell(f"cell_A{i}", loc)

locsB = [
    (1, 0, 0),
    (3, 0, 0),
    (1, 2, 0),
    (3, 2, 0),
    (1, 0, -2),
    (3, 0, -2),
    (1, 2, -2),
    (3, 2, -2),
]

for i, loc in enumerate(locsB):
    cellsB.create_cell(f"cell_B{i}", loc, color=(0.1, 0, 0))

sim = Simulator([cellsA, cellsB])
sim.setup_world()
sim.add_handlers(
    [
        GrowthPIDHandler(target_volume=30),
        RecenterHandler(),
        SceneExtensionHandler(end=300),
    ]
)

sim.render(end=300)
