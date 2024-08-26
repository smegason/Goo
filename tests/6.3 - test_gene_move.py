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
cell = celltype.create_cell("cell", (1, -1, 0))
cell.pressure = 1
# cell.move((0, 0, -1))

cell.link_gene_to_property(mol, "motion_direction")

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

points = [[0, 0, 1], [0, 0, -1], [0, 1, 0], [0, -1, 0], [1, 0, 0], [-1, 0, 0]]
for point in points:
    print("Point", point, "\tConcentration:", diffsys.get_concentration(mol, point))

for i in range(1, 11):
    bpy.context.scene.frame_set(i)
    print("Motion direction:", cell.motion_force.loc - cell.loc)
    print("New concentration", cell.metabolites[mol])

# sim.run(1)
# sim.render_animation(end=60)
