from importlib import reload
import goo
from goo.handler import *
from goo.utils import PhysicsConstructor, ClothConstructor

reload(goo)
goo.reset_modules()
goo.reset_scene()

celltype = goo.SimpleType("cellsA")
for i in range(10):
    cell = celltype.create_cell(
        f"cell{i}",
        (0, 0, 0),
        size=2,
        physics_constructor=PhysicsConstructor(
            ClothConstructor,
        ),
    )
    cell.stiffness = 15
    cell.pressure = 1

sim = goo.Simulator([celltype])
sim.setup_world(seed=i)
sim.add_handlers(
    [
        GrowthPIDHandler(),
        AdhesionLocationHandler(),
        RandomMotionHandler(ForceDist.CONSTANT, max_strength=8000),
        RemeshHandler(),
        DataExporter(
            path=f"/tmp/motion_out{i}.json",
            options=DataFlag.TIMES | DataFlag.MOTION_PATH | DataFlag.FORCE_PATH,
        ),
    ]
)

sim.run(end=100)
