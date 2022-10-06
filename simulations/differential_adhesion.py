from Goo import Goo
import bpy
import numpy as np
import importlib

importlib.reload(Goo)

Goo.setup_world()

cells = []
forces = []

c1 = Goo.Cell("c1", loc=(0, 0, 0))
Goo.make_cell(c1)
cells.append(c1)

c2 = Goo.Cell("c2", loc=(2, 2, 2))
Goo.make_cell(c2)
cells.append(c2)

c3 = Goo.Cell("c3", loc=(-2, -2, -2))
Goo.make_cell(c3)
cells.append(c3)

c4 = Goo.Cell("c4", loc=(-4, -4, -4))
Goo.make_cell(c4)
cells.append(c4)


f1 = Goo.Force("f1", "c1", -800)
Goo.make_force(f1)
forces.append(f1)

f2 = Goo.Force("f2", "c2", -800)
Goo.make_force(f2)
forces.append(f2)

f3 = Goo.Force("f3", "c1", -200)
Goo.make_force(f3)
forces.append(f3)

f4 = Goo.Force("f4", "c2", -200)
Goo.make_force(f4)
forces.append(f4)

f5 = Goo.Force("f5", "c3", -800)
Goo.make_force(f5)
forces.append(f5)

f6 = Goo.Force("f6", "c4", -800)
Goo.make_force(f6)
forces.append(f6)

f7 = Goo.Force("f7", "c3", -200)
Goo.make_force(f7)
forces.append(f7)

f8 = Goo.Force("f8", "c4", -200)
Goo.make_force(f8)
forces.append(f8)


def force_handler(scene):
    for i in range(len(forces)):
        assoc_cell = forces[i].associated_cell
        bpy.context.view_layer.objects.active = bpy.data.objects[assoc_cell]
        dg = bpy.context.evaluated_depsgraph_get()
        cell_eval = bpy.data.objects[assoc_cell].evaluated_get(dg)
        vertices = cell_eval.data.vertices
        vert_coords = [(cell_eval.matrix_world @ v.co) for v in vertices]
        vert_coords = np.asarray(vert_coords)

        x = vert_coords[:, 0]
        y = vert_coords[:, 1]
        z = vert_coords[:, 2]
        COM = (np.mean(x), np.mean(y), np.mean(z))
        bpy.data.objects[forces[i].name].location = COM


bpy.app.handlers.frame_change_post.clear()
bpy.app.handlers.frame_change_post.append(force_handler)
