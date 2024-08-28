from importlib import reload
import goo

reload(goo)
goo.reset_modules()
goo.reset_scene()

# Defining cells
celltype = goo.CellType("cellA")
celltype.homo_adhesion_strength = 1000
cell = celltype.create_cell(name="cell", loc=(0, 0, 0), color=(0, 1, 1))
cell.stiffness = 2
cell.pressure = 5

sim = goo.Simulator(celltypes=[celltype], time=500, physics_dt=1)
sim.setup_world()
sim.add_handlers(
    [
        goo.GrowthPIDHandler(target_volume=70),
        goo.SizeDivisionHandler(goo.BisectDivisionLogic, mu=60, sigma=2),
        goo.AdhesionLocationHandler(),
        goo.RemeshHandler(),
    ]
)