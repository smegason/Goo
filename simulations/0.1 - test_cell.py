# Goal: test the basic methods of calculating cell properties,
# comparing them to the original implementations.

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

print("Volume", cell.volume())
print("COM:", cell.COM())

print("---- Long axis ----")
print(g.get_long_axis(cell.obj))
axis = cell.major_axis()
print(axis.axis(), axis.length(), axis.endpoints(local_coords=True))

print("---- Short axis ----")
print(g.get_minor_axis(cell.obj))
m_axis = cell.minor_axis()
print(
    m_axis.axis(),
    m_axis.length(local_coords=True),
    m_axis.endpoints(),
)
