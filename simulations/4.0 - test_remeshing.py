from importlib import reload
import goo
from goo.division import *
from goo.handler import *

from scipy.spatial.distance import pdist, cdist, squareform

reload(goo)
goo.reset_modules()
goo.reset_scene()

celltype = goo.CellType("A", physics_enabled=False)
cell = celltype.create_cell("cellA", (0, 0, 0))

np.random.seed(1)
for i in range(10):
    loc = np.random.random_sample(3) * 6 - 3
    celltype.create_cell(f"cellA{i}", loc)

# model_mesh = cell.obj.data.copy()
# model_obj = bpy.data.objects.new("model", model_mesh)
# bpy.context.scene.collection.objects.link(model_obj)
# mod = model_obj.modifiers.new("Shrinkwrap", "SHRINKWRAP")
# mod.wrap_method = "NEAREST_SURFACEPOINT"
# mod.wrap_mode = "ON_SURFACE"

# logic = BisectDivisionLogic()
# mother, daughter = cell.divide(logic)
# logic.flush()
# mother.remesh(0.1)
# daughter.remesh(0.1)

# mother.obj.hide_set(True)
# daughter.obj.hide_set(True)

# mod.target = mother.obj
# dg = bpy.context.evaluated_depsgraph_get()
# mother.data = model_obj.evaluated_get(dg).copy()


def _contact_area(cell1: Cell, cell2: Cell, threshold=0.1):
    fs1 = cell1.obj.data.polygons
    fs2 = cell2.obj.data.polygons

    centers1 = [cell1.obj.matrix_world @ f.center for f in fs1]
    centers2 = [cell2.obj.matrix_world @ f.center for f in fs2]

    dists = np.array(cdist(centers1, centers2, "euclidean"))

    contact_fs1 = np.any(dists < threshold, axis=1)
    contact_fs2 = np.any(dists < threshold, axis=0)

    areas1 = np.array([f.area for f in fs1])
    areas2 = np.array([f.area for f in fs2])

    ratio1 = np.sum(areas1[contact_fs1]) / np.sum(areas1)
    ratio2 = np.sum(areas2[contact_fs2]) / np.sum(areas2)

    return ratio1, ratio2


def _contact_areas(cells, threshold):
    coms = [cell.COM() for cell in cells]
    dists = squareform(pdist(coms, "euclidean"))
    print(dists)

    mask = dists < threshold
    mask = np.triu(mask, k=1)

    pairs = np.where(mask)

    areas = {cell.name: [] for cell in cells}
    for i, j in zip(pairs[0], pairs[1]):
        ratio_i, ratio_j = _contact_area(cells[i], cells[j])
        areas[cells[i].name].append((cells[j].name, ratio_i))
        areas[cells[j].name].append((cells[i].name, ratio_j))

    return areas


out = _contact_areas(celltype.cells, 2)
print(out)
