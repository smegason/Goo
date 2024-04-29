from importlib import reload
import goo
from goo import *

reload(goo)
reset_modules()
reset_scene()

cellsA = SimpleType("A")
cellsA.homo_adhesion_strength = 2500

cellsA.create_cell("A1", (0, 0, 0))

sim = Simulator([cellsA])
sim.setup_world()
sim.add_handlers(
    [
        GrowthPIDHandler(target_volume=30),
        AdhesionLocationHandler(),
        TimeDivisionHandler(BisectDivisionLogic, mu=30),
        RemeshHandler(),
        DataExporter(
            path="/tmp/out.json",
            options=DataFlag.TIMES | DataFlag.VOLUMES | DataFlag.PRESSURES,
        ),
    ]
)
sim.run(179)
