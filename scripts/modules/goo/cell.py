from functools import reduce
from typing import Optional, Union, Callable
import numpy as np

import bpy, bmesh
from bpy.types import Modifier, ClothModifier, CollisionModifier
from mathutils import *

from goo.force import *
from goo.utils import *


class Cell(BlenderObject):
    def __init__(self, obj: bpy.types.Object, mat=None):
        super(Cell, self).__init__(obj)

        # Set up effector collections
        self._effectors = bpy.data.collections.new(f"{obj.name}_effectors")
        bpy.context.scene.collection.children.link(self._effectors)
        self.add_effector(ForceCollection.global_forces())  # link global forces

        self._mat = mat
        self.obj.data.materials.append(mat)

        self.celltype: CellType = None

        self._physics_enabled = False
        self.mod_settings = []

        self._homo_adhesion: AdhesionForce = None
        self._hetero_adhesions: list[AdhesionForce] = []
        self._motion_force: MotionForce = None

    @property
    def name(self) -> str:
        return self.obj.name

    @name.setter
    def name(self, name: str):
        old_name = self.obj.name

        self.obj.name = name
        self.obj.data.name = f"{name}_mesh"
        self._effectors.name = f"{name}_effectors"
        if self._mat:
            self._mat.name = f"{name}_material"

        for force in self._hetero_adhesions:
            force.name = force.name.replace(old_name, name, 1)
        if self.motion_force:
            self.motion_force.name = self.motion_force.name.replace(old_name, name, 1)

    def copy(self) -> "Cell":
        """Copies the cell's underlying Blender object, object data, and material if applicable.
        Returns a new cell with the copied data.

        Note that any custom collections will not be updated. It is advised that
        physics modifiers are set up again for the daughter cell after calling `copy`.
        """
        # Set up object, mesh, and material for copy
        obj_copy = self.obj.copy()
        obj_copy.data = self.obj.data.copy()
        bpy.context.scene.collection.objects.link(obj_copy)

        if self._mat is not None:
            obj_copy.data.materials.clear()
            mat_copy = self._mat.copy()
        else:
            mat_copy = None

        cell_copy = Cell(obj_copy, mat_copy)
        cell_copy._physics_enabled = self.physics_enabled
        cell_copy._update_cloth()

        return cell_copy

    # ----- CUSTOM PROPERTIES -----
    def __setitem__(self, k: str, v: Union[float, list[float], int, list[int], str]):
        self.obj.data[k] = v

    def __contains__(self, k: str):
        return k in self.obj.data.keys()

    def __getitem__(self, k):
        return self.obj.data[k]

    # ----- BASIC FUNCTIONS -----
    @property
    def obj_eval(self):
        dg = bpy.context.evaluated_depsgraph_get()
        obj_eval = self.obj.evaluated_get(dg)
        return obj_eval

    def vertices(self, local_coords=False):
        verts = self.obj_eval.data.vertices
        if local_coords:
            return [v.co for v in verts]
        else:
            return [self.obj_eval.matrix_world @ v.co for v in verts]

    def volume(self) -> float:
        """Calculates the volume of the Blender mesh."""
        bm = bmesh.new()
        bm.from_mesh(self.obj_eval.to_mesh())
        bm.transform(self.obj_eval.matrix_world)
        volume = bm.calc_volume()
        bm.free()

        return volume

    def COM(self, local_coords: bool = False) -> Vector:
        """Calculates the center of mass of a mesh."""
        vert_coords = self.vertices(local_coords)
        com = Vector(np.mean(vert_coords, axis=0))
        return com

    def _get_eigenvector(self, n):
        """Returns the nth eigenvector (axis) in object space as a line defined by two Vectors."""
        verts = self.vertices()

        # Calculate the eigenvectors and eigenvalues of the covariance matrix
        covariance_matrix = np.cov(verts, rowvar=False)
        eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)
        eigenvectors = eigenvectors[:, eigenvalues.argsort()[::-1]]
        axis = eigenvectors[:, n]

        # sort indices by distance
        vert_indices = np.argsort(np.asarray(verts).dot(axis))
        first_vertex = verts[vert_indices[0]]
        last_vertex = verts[vert_indices[-1]]

        return Axis(Vector(axis), first_vertex, last_vertex, self.obj_eval.matrix_world)

    def major_axis(self) -> Axis:
        return self._get_eigenvector(0)

    def minor_axis(self) -> Axis:
        return self._get_eigenvector(1)

    def divide(self, divisionLogic) -> tuple["Cell", "Cell"]:
        mother, daughter = divisionLogic.make_divide(self)
        if mother.celltype:
            mother.celltype.add_cell(daughter)
        return mother, daughter

    def recenter(self):
        """Recenter origin to COM."""
        bm = bmesh.new()
        bm.from_mesh(self.obj_eval.to_mesh())

        com = self.COM()
        bmesh.ops.translate(bm, verts=bm.verts, vec=-self.COM(local_coords=True))

        bm.to_mesh(self.obj.data)
        bm.free()

        self.loc = com

    def remesh(self, voxel_size=0.25, smooth: bool = True):
        # use of object ops is 2x faster than remeshing with modifiers
        self.obj.data.remesh_mode = "VOXEL"
        self.obj.data.remesh_voxel_size = voxel_size
        with bpy.context.temp_override(active_object=self.obj, object=self.obj):
            bpy.ops.object.voxel_remesh()

        for f in self.obj.data.polygons:
            f.use_smooth = smooth

    def recolor(self, color):
        r, g, b = color
        _, _, _, a = self._mat.diffuse_color
        self._mat.diffuse_color = (r, g, b, a)

        if self._mat.use_nodes:
            for node in self._mat.node_tree.nodes:
                if "Base Color" in node.inputs:
                    _, _, _, a = node.inputs["Base Color"].default_value
                    node.inputs["Base Color"].default_value = r, g, b, a

    # ----- PHYSICS -----
    def get_modifier(self, type) -> Optional[Modifier]:
        return next((m for m in self.obj.modifiers if m.type == type), None)

    @property
    def cloth_mod(self) -> Optional[ClothModifier]:
        return self.get_modifier("CLOTH")

    @property
    def collision_mod(self) -> Optional[CollisionModifier]:
        return self.get_modifier("COLLISION")

    @property
    def physics_enabled(self) -> bool:
        return self._physics_enabled

    def _update_cloth(self):
        """Ensure that cloth effectors setting is correctly set."""
        if self.cloth_mod:
            self.cloth_mod.settings.effector_weights.collection = self._effectors

    def setup_physics(self, physics_constructor: PhysicsConstructor):
        physics_constructor(self.obj)
        self._update_cloth()
        self._physics_enabled = True

    def enable_physics(self):
        if self._physics_enabled:
            raise RuntimeError(
                f"{self.name}: physics must be disabled before enabling."
            )
        # recreate modifier stack
        for name, type, settings in self.mod_settings:
            mod = self.obj.modifiers.new(name=name, type=type)
            declare_settings(mod, settings)
        self.mod_settings.clear()

        # ensure cloth mod is set correctly
        self._update_cloth()

        for force in self.adhesion_forces:
            force.enable()
        self._physics_enabled = True

    def disable_physics(self):
        if not self._physics_enabled:
            raise RuntimeError(
                f"{self.name}: physics must be set up and/or enabled before disabling."
            )
        for mod in self.obj.modifiers:
            name, type = mod.name, mod.type
            settings = store_settings(mod)
            self.mod_settings.append((name, type, settings))
            self.obj.modifiers.remove(mod)

        for force in self.adhesion_forces:
            force.disable()
        self._physics_enabled = False

    @property
    def stiffness(self) -> float:
        return self.cloth_mod.settings.tension_stiffness

    @stiffness.setter
    def stiffness(self, stiffness: float):
        self.cloth_mod.settings.tension_stiffness = stiffness
        self.cloth_mod.settings.compression_stiffness = stiffness
        self.cloth_mod.settings.shear_stiffness = stiffness

    @property
    def pressure(self) -> float:
        return self.cloth_mod.settings.uniform_pressure_force

    @pressure.setter
    def pressure(self, pressure: float):
        self.cloth_mod.settings.uniform_pressure_force = pressure

    # ----- FORCES -----
    def add_effector(self, force: Force | ForceCollection):
        """Add a force or a collection of forces that affects this cell."""
        if isinstance(force, Force):
            self._effectors.objects.link(force.obj)
        elif isinstance(force, ForceCollection):
            self._effectors.children.link(force.collection)

    def remove_effector(self, force: Force | ForceCollection):
        """Removes a force or a collection of forces that affects this cell."""
        if isinstance(force, Force):
            self._effectors.objects.unlink(force.obj)
        elif isinstance(force, ForceCollection):
            self._effectors.children.unlink(force.collection)

    def link_adhesion_force(self, force: AdhesionForce):
        """Add force that stems from this cell."""
        self._hetero_adhesions.append(force)

    @property
    def homo_adhesion(self) -> AdhesionForce:
        return self._homo_adhesion

    @homo_adhesion.setter
    def homo_adhesion(self, force: AdhesionForce):
        self._homo_adhesion = force

    @property
    def adhesion_forces(self) -> list[AdhesionForce]:
        return [self._homo_adhesion] + self._hetero_adhesions

    @property
    def motion_force(self) -> Force:
        return self._motion_force

    @motion_force.setter
    def motion_force(self, force: MotionForce):
        if self.motion_force:
            self.remove_effector(self.motion_force)
        self.add_effector(force)
        self._motion_force = force

    def move_towards(self, dir: Vector):
        """Sets motion force a certain direction."""
        if not self._motion_force:
            raise RuntimeError(
                f"Cell {self.name} does not have an associated motion force!"
            )
        motion_loc = self.loc + dir.normalized() * (2 + self.major_axis().length())
        self._motion_force.set_loc(motion_loc, self.loc)


def create_cell(
    name, loc, color=None, physics_constructor=None, physics_on=True, **kwargs
) -> Cell:
    cell = CellType.default_celltype().create_cell(
        name, loc, color, physics_constructor, physics_on, **kwargs
    )
    return cell


def store_settings(mod: bpy.types.bpy_struct):
    settings = {}
    for p in mod.bl_rna.properties:
        id = p.identifier
        if not p.is_readonly:
            settings[id] = getattr(mod, id)
        elif id in ["settings", "collision_settings", "effector_weights"]:
            settings[id] = store_settings(getattr(mod, id))
    return settings


def declare_settings(mod: bpy.types.bpy_struct, settings: dict):
    for id, setting in settings.items():
        if isinstance(setting, dict):
            declare_settings(getattr(mod, id), settings[id])
        else:
            setattr(mod, id, setting)


class CellType:
    mesh_kwargs = {}
    physics_constructor = PhysicsConstructor(
        SubsurfConstructor,
        ClothConstructor,
        CollisionConstructor,
        RemeshConstructor,
    )
    color = (0.007, 0.021, 0.3)
    _default_celltype = None

    def __init__(self, name, physics_enabled=True):
        self._homo_adhesions = ForceCollection(name)
        self._cells = set()

        self._physics_enabled = physics_enabled
        self._hetero_adhesions = {}

        self.homo_adhesion_strength = 2000
        self.motion_strength = 0

    @staticmethod
    def default_celltype():
        if CellType._default_celltype is None:
            CellType._default_celltype = CellType("default")
        return CellType._default_celltype

    @property
    def name(self):
        return self._homo_adhesions.name

    def add_cell(self, cell: Cell):
        """Adds cell to the cell type, activating its physics and constructing
        appropriate adhesion and motion forces."""
        cell.celltype = self
        self._cells.add(cell)

        if not self._physics_enabled:
            return

        # add homotypic adhesion force
        homo_adhesion = get_adhesion(self.homo_adhesion_strength, obj=cell.obj)
        cell.homo_adhesion = homo_adhesion

        self._homo_adhesions.add_force(homo_adhesion)
        cell.add_effector(self._homo_adhesions)

        # add heterotypic adhesion forces
        for celltype, item in self._hetero_adhesions.items():
            outgoing_forces, incoming_forces, strength = item
            hetero_adhesion = get_adhesion(
                strength, name=f"{cell.name}_to_{celltype.name}", loc=cell.loc
            )
            cell.add_effector(incoming_forces)

            outgoing_forces.add_force(hetero_adhesion)
            cell.link_adhesion_force(hetero_adhesion)

        # add motion force
        motion = get_motion(
            name=f"{cell.name}_motion",
            loc=(0, 0, 0),
            strength=self.motion_strength,
        )
        cell.motion_force = motion

    # TODO: implement correctly
    def _remove_cell(self, cell: Cell):
        pass
        # cell.celltype = None
        # self._cells.remove(cell)

    def create_cell(
        self,
        name,
        loc,
        color=None,
        physics_constructor=None,
        physics_on=True,
        **mesh_kwargs,
    ) -> Cell:
        """Creates Blender object and materials, and creates a cell.

        Note: `create_cell` does not directly activate physics, which is handled by the
        `add_cell` function. This is significant for division, which will use the default
        CellType settings to set up cell physics, rather than any custom settings of the
        dividing cell."""
        mesh_kwargs = dict(self.__class__.mesh_kwargs, **mesh_kwargs)
        if color is None:
            color = self.__class__.color

        obj = create_mesh(name, loc, mesh="icosphere", **mesh_kwargs)
        bpy.context.scene.collection.objects.link(obj)

        mat = create_material(f"{name}_material", color=color) if color else None
        cell = Cell(obj, mat)
        cell.remesh()

        # enable physics for cell
        if physics_constructor is None:
            physics_constructor = self.__class__.physics_constructor
        if self._physics_enabled and physics_on:
            cell.setup_physics(physics_constructor)

        self.add_cell(cell)
        return cell

    # TODO: create cells in shape or in specified locations of specified cell types
    def create_cells(celltype, locs=None, shape=None, size=None):
        pass

    @property
    def cells(self) -> list[Cell]:
        return list(self._cells)

    def set_hetero_adhesion(self, other_celltype: "CellType", strength: float):
        outgoing_forces = ForceCollection(f"{self.name}_to_{other_celltype.name}")
        incoming_forces = ForceCollection(f"{other_celltype.name}_to_{self.name}")
        self._hetero_adhesions[other_celltype] = (
            outgoing_forces,
            incoming_forces,
            strength,
        )
        other_celltype._hetero_adhesions[self] = (
            incoming_forces,
            outgoing_forces,
            strength,
        )


class SimpleType(CellType):
    color = None
    physics_constructor = PhysicsConstructor(
        ClothConstructor,
        CollisionConstructor,
    )


class YolkType(CellType):
    mesh_kwargs = {
        "size": 10,
        "subdivisions": 4,
    }
    physics_constructor = PhysicsConstructor(
        YolkClothConstructor,
        CollisionConstructor,
    )
    color = (0.64, 0.64, 0.64)
