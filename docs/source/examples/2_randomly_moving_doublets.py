from importlib import reload
import goo

reload(goo)
goo.reset_modules()
goo.reset_scene()

cellsA = goo.CellType("A")
cellsB = goo.CellType("B")
cellsC = goo.CellType("C")

cellsA.homo_adhesion_strength = 15
cellsA.stiffness = 3  # ratio = 5
cellsB.homo_adhesion_strength = 100
cellsB.stiffness = 2  # ratio = 50
cellsC.homo_adhesion_strength = 500
cellsC.stiffness = 1  # ratio = 500

cellsA.create_cell("A1", (-10, +1.75, 0), color=(1, 1, 0), size=1.6)
cellsA.create_cell("A2", (-10, -1.75, 0), color=(1, 1, 0), size=1.6)

cellsB.create_cell("B1", (0, +1.75, 0), color=(0, 1, 1), size=1.6)
cellsB.create_cell("B2", (0, -1.75, 0), color=(0, 1, 1), size=1.6)

cellsC.create_cell("C1", (10, +1.75, 0), color=(1, 0, 1), size=1.6)
cellsC.create_cell("C2", (10, -1.75, 0), color=(1, 0, 1), size=1.6)

sim = goo.Simulator([cellsA, cellsB, cellsC], time=130, physics_dt=1)
sim.setup_world(seed=2024)
sim.add_handlers(
    [
        goo.GrowthPIDHandler(target_volume=50),  # in um3
        goo.RecenterHandler(),  # no parameters needed
        goo.SizeDivisionHandler(goo.BisectDivisionLogic, mu=50, sigma=1),  # in um3
        goo.RandomMotionHandler(goo.ForceDist.GAUSSIAN, max_strength=1000),
    ]
)
