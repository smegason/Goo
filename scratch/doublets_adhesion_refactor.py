from importlib import reload
import goo
from goo import *
from goo.utils import *
from math import pi

reload(goo)
reset_modules()
reset_scene()


class CastConstructor(ModConstructor):
    name = "Cast"
    type = "CAST"

    def setup_mod(self, mod: bpy.types.CastModifier):
        mod.cast_type = "SPHERE"
        mod.factor = 0.25
        mod.radius = 0
        mod.size = 1.65
        mod.use_radius_as_size = False


class CastType(SimpleType):
    physics_constructor = PhysicsConstructor(
        CastConstructor,
        ClothConstructor,
        CollisionConstructor,
    )


cellsA = CastType("A")
cellsB = CastType("B")
cellsC = CastType("C")

cellsA.homo_adhesion_strength = 5000
cellsB.homo_adhesion_strength = 7500
cellsC.homo_adhesion_strength = 10000

# cellsA.set_hetero_adhesion(cellsB, 10000)
# cellsA.set_hetero_adhesion(cellsC, 10000)


cellsA.create_cell("A1", (1.75, 0, 0), color=(0, 0, 0.1), size=1.6)
cellsA.create_cell("A2", (-1.75, 0, 0), color=(0, 0, 0.1), size=1.6)

cellsB.create_cell("B1", (1.75, 3.2, 0), color=(0, 0.1, 0.6), size=1.6)
cellsB.create_cell("B2", (-1.75, 3.2, 0), color=(0, 0.1, 0.6), size=1.6)

cell1 = cellsC.create_cell("C1", (1.75, 6.4, 0), size=1.6)
cell2 = cellsC.create_cell("C2", (-1.75, 6.4, 0), size=1.6)

sim = Simulator([cellsA, cellsB, cellsC])
sim.setup_world()
sim.add_handlers(
    [
        GrowthPIDHandler(target_volume=20),
        RecenterHandler(),
        # RemeshHandler(),
        SceneExtensionHandler(end=500),
    ]
)

# for i in range(50):
#     bpy.context.scene.frame_set(i + 1)
#     print(cell1.COM(), cell2.COM())
