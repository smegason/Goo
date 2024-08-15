from importlib import reload
import goo
from goo.cell import *
from goo.handler import *
from goo.simulator import Render

reload(goo)
goo.reset_modules()
goo.reset_scene()

mol = Molecule("mol", conc=0.01, D=1, gradient="linear")
diffsys = DiffusionSystem(molecules=[mol])

celltype = CellType("cellsA", pattern="simple")
celltype.motion_strength = 50
cell = celltype.create_cell("cell", (0, 0, 0))
cell.pressure = 1
cell.move((0, 0, -1))

cell.link_property_to_gene("motion_force", mol, pscale=(0, 50), gscale=(1, 1.5))

sim = goo.Simulator(cells=[cell], diffsystem=diffsys)
sim.setup_world()
sim.add_handlers(
    [
        RecenterHandler(),
        DiffusionHandler(),
        NetworkHandler(),
        ColorizeHandler(Colorizer.GENE, mol, range=(1, 1.5)),
    ]
)

sim.render_animation(end=60)
