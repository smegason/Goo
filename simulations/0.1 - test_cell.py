# Goal: test the creation of a world, and population of cells within that world
# through base Cell creation and CellType creation.

from importlib import reload
import goo
from goo import goo as g

reload(goo)
goo.reset_modules()
goo.reset_scene()

cell = goo.create_cell(
    "cell",
    (1, 2, 3),
    rotation=(1, 2, 3),
    scale=(1, 2, 3),
    subdivisions=4,
    physics_on=False,
)
print(cell.obj.location)

print("Volume", cell.get_volume())
print("COM:", cell.get_COM())

print("---- Long axis ----")
print(g.get_long_axis(cell.obj))
axis = cell.get_major_axis()
print(axis.axis(), axis.length(), axis.endpoints(local_coords=True))

print("---- Short axis ----")
print(g.get_minor_axis(cell.obj))
m_axis = cell.get_minor_axis()
print(
    m_axis.axis(),
    m_axis.length(local_coords=True),
    m_axis.endpoints(),
)

plane = cell.create_division_plane()
plane.hide_set(False)
