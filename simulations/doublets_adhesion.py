from importlib import reload
import goo
from goo import *

reload(goo)
reset_modules()
reset_scene()

cellsA = CellType("A")
cellsB = CellType("B")
cellsC = CellType("C")

cellsA.homo_adhesion_strength = 2500
cellsB.homo_adhesion_strength = 5000
cellsC.homo_adhesion_strength = 7500

# Uncomment for heteroadhesions
# cellsA.set_hetero_adhesion(cellsB, 10000)
# cellsA.set_hetero_adhesion(cellsC, 10000)
# cellsB.set_hetero_adhesion(cellsC, 10000)

cellsA.create_cell("A1", (+1.75, -5, 0), color=(0.5, 0, 0), size=1.6)
cellsA.create_cell("A2", (-1.75, -5, 0), color=(0.5, 0, 0), size=1.6)

# cellsB.create_cell("B1", (+1.75, 0, 0), color=(0, 0.5, 0), size=1.6)
# cellsB.create_cell("B2", (-1.75, 0, 0), color=(0, 0.5, 0), size=1.6)

# cellsC.create_cell("C1", (+1.75, 5, 0), color=(0, 0, 0.5), size=1.6)
# cellsC.create_cell("C2", (-1.75, 5, 0), color=(0, 0, 0.5), size=1.6)

sim = Simulator([cellsA, cellsB, cellsC])
sim.setup_world()
sim.add_handlers(
    [
        GrowthPIDHandler(target_volume=20),
        AdhesionLocationHandler(),
        # RemeshHandler(),
        DataExporter(
            path="/tmp/out.json", options=DataFlag.TIMES | DataFlag.CONTACT_AREAS
        ),
    ]
)
# sim.run(20)
