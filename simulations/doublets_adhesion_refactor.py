from importlib import reload
import goo
from goo import *

reload(goo)
reset_modules()
reset_scene()

cellsA = SimpleType("A")
cellsB = SimpleType("B")
cellsC = SimpleType("C")

cellsA.homo_adhesion_strength = 5000
cellsB.homo_adhesion_strength = 7500
cellsC.homo_adhesion_strength = 10000

cellsA.set_hetero_adhesion(cellsB, 10000)
cellsA.set_hetero_adhesion(cellsC, 10000)

cellsA.create_cell("A1", (1.75, 0, 0), color=(0, 0, 0.1), scale=(1.6, 1.6, 1.6))
cellsA.create_cell("A2", (-1.75, 0, 0), color=(0, 0, 0.1), scale=(1.6, 1.6, 1.6))

cellsB.create_cell("B1", (1.75, 3.2, 0), color=(0, 0.1, 0), scale=(1.6, 1.6, 1.6))
cellsB.create_cell("B2", (-1.75, 3.2, 0), color=(0, 0.1, 0), scale=(1.6, 1.6, 1.6))

cellsC.create_cell("C1", (1.75, 6.4, 0), color=(0.1, 0, 0), scale=(1.6, 1.6, 1.6))
cellsC.create_cell("C2", (-1.75, 6.4, 0), color=(0.1, 0, 0), scale=(1.6, 1.6, 1.6))

sim = Simulator([cellsA, cellsB, cellsC])
sim.setup_world()
sim.add_handlers(
    [
        GrowthPIDHandler(target_volume=20),
        AdhesionLocationHandler(),
        SceneExtensionHandler(end=500),
    ]
)
