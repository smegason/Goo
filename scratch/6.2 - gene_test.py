from importlib import reload
import goo
from goo.cell import *
from goo.gene import *
from goo.division import *
from goo.handler import *

reload(goo)
goo.reset_modules()
goo.reset_scene()

# Defining cells
x = Gene("x")
y = Gene("y")
z = Gene("z")


class Pattern(SimplePattern):
    circuits = [
        DegFirstOrder(x, 0.1),
        DegFirstOrder(y, 0.1),
        DegFirstOrder(z, 0.1),
        ProdRepression(y, x, kcat=0.4, n=3),
        ProdRepression(z, y, kcat=0.4, n=3),
        ProdRepression(x, z, kcat=0.4, n=3),
    ]


celltype = CellType("cellA", pattern=Pattern())

cell1 = celltype.create_cell(
    name="cell1",
    loc=(0, -1.5, 0),
    physics_enabled=False,
    growth_enabled=False,
)
cell1.metabolites = {x: 2, y: 0.1, z: 0.1}

cell2 = celltype.create_cell(
    name="cell2",
    loc=(0, 1.5, 0),
    physics_enabled=False,
    growth_enabled=False,
)
cell2.metabolites = {x: 0.1, y: 2, z: 0.1}

cell3 = celltype.create_cell(
    name="cell2",
    loc=(0, 4.5, 0),
    physics_enabled=False,
    growth_enabled=False,
)
cell3.metabolites = {x: 0.1, y: 0.1, z: 2}


sim = goo.Simulator(celltypes=[celltype], time=250)
sim.setup_world()
sim.add_handlers([NetworkHandler(), ColorizeHandler(Colorizer.GENE, x, range=(0, 2))])
# sim.render()

print("Done loading")
