from importlib import reload
import goo


reload(goo)
goo.reset_modules()
goo.reset_scene()

# Defining cells
celltype = goo.CellType("A", target_volume=70, pattern="simple")
celltype.homo_adhesion_strength = 500

cell1 = celltype.create_cell(name="cell1", loc=(0, 0, 0), color=(0, 1, 1))
cell1.stiffness = 2
cell1.pressure = 5

sim = goo.Simulator([celltype], time=200, physics_dt=1)
sim.setup_world(seed=2025)
sim.add_handlers(
    [
        goo.GrowthPIDHandler(),                                             # in um3
        goo.RecenterHandler(),
        goo.RemeshHandler(),
        goo.SizeDivisionHandler(goo.BisectDivisionLogic, mu=65, sigma=2),   # in um3
    ]
)
