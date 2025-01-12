from importlib import reload
import goo


reload(goo)
goo.reset_modules()
goo.reset_scene()

# Defining cells
celltype = goo.SimpleType("cellA")
celltype.homo_adhesion_strength = 1000
cell1 = celltype.create_cell(name="cell", loc=(0, 0, 0), color=(0, 1, 1))
cell1.stiffness = 2
cell1.pressure = 5
output_dir = "<output/directory/>"

sim = goo.Simulator(celltypes=[celltype], time=500, physics_dt=1)
sim.setup_world()
sim.add_handlers(
    [
        goo.GrowthPIDHandler(target_volume=50),
        goo.SizeDivisionHandler(goo.BisectDivisionLogic, mu=50, sigma=2),
        goo.RecenterHandler(),
        goo.RemeshHandler(),
        goo.SliceExporter(output_dir=output_dir, z_range=(3, -3), z_step=0.2),
    ]
)
