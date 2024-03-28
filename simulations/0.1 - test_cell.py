# Goal: test the creation of a world, and population of cells within that world
# through base Cell creation and CellType creation.

from importlib import reload
import goo
from goo import goo as g
from goo.cell import *
from goo.reloader import *

reload(goo)
reset_modules()
reset_scene()

cell = create_cell("cell", (1, 2, 3), rotation=(1, 2, 3), scale=(2, 1, 0.5))

print("Volume", cell.volume())
print("COM:", cell.COM())

print("---- Long axis ----")
print(g.get_long_axis(cell.obj))
axis = cell.major_axis()
print(axis.axis(), axis.length(global_coords=True), axis.endpoints())

print("---- Short axis ----")
print(g.get_minor_axis(cell.obj))
m_axis = cell.minor_axis()
print(m_axis.axis(global_coords=True), m_axis.length(), m_axis.endpoints())

plane = cell.create_division_plane()
plane.hide_set(False)
