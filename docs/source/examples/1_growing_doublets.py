from importlib import reload
import goo

reload(goo)
goo.reset_modules()
goo.reset_scene()

cellsA = goo.CellType("A", target_volume=60, pattern="simple")
cellsB = goo.CellType("B", target_volume=60, pattern="simple")
cellsC = goo.CellType("C", target_volume=60, pattern="simple")

cellsA.homo_adhesion_strength = 1
cellsA.stiffness = 1  # ratio = 1
cellsB.homo_adhesion_strength = 500
cellsB.stiffness = 1  # ratio = 500
cellsC.homo_adhesion_strength = 2000
cellsC.stiffness = 1  # ratio = 2000

cellsA.create_cell("A1", (-20, +1.75, 0), color=(1, 1, 0), size=1.6)
cellsA.create_cell("A2", (-20, -1.75, 0), color=(1, 1, 0), size=1.6)

cellsB.create_cell("B1", (0, +1.75, 0), color=(0, 1, 1), size=1.6)
cellsB.create_cell("B2", (0, -1.75, 0), color=(0, 1, 1), size=1.6)

cellsC.create_cell("C1", (20, +1.75, 0), color=(1, 0, 1), size=1.6)
cellsC.create_cell("C2", (20, -1.75, 0), color=(1, 0, 1), size=1.6)

sim = goo.Simulator([cellsA, cellsB, cellsC], time=180, physics_dt=1)
sim.setup_world(seed=2024)
sim.add_handlers(
    [
        goo.GrowthPIDHandler(),                                             # in um3
        goo.RecenterHandler(),
        goo.SizeDivisionHandler(goo.BisectDivisionLogic, mu=50, sigma=1),   # in um3
        goo.RandomMotionHandler(goo.ForceDist.GAUSSIAN, strength=500),
    ]
)
