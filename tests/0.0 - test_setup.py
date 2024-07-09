# Goal: test the ability to import and reload modules

from importlib import reload
import goo

reload(goo)
goo.reset_modules()
goo.reset_scene()

cell = goo.create_cell("cell", (0, 0, 0), physics_on=False)
