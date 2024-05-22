from functools import reduce
from typing import Optional, Union, Callable
import numpy as np

import bpy, bmesh
from bpy.types import Modifier, ClothModifier, CollisionModifier
from mathutils import *

from goo.force import *
from goo.utils import *


class Cell(BlenderObject):
    """A cell.

    Cells are represented in Blender by mesh objects. They can interact with
    physics by adding Blender modifiers, and the forces that influence its
    motion is determined by an associated collection.

    Args:
        obj: Blender object to be used as the representation of the cell.

    Attributes:
        celltype (CellType): The cell type to which the cell belongs.
    """

    def __init__(self, obj: bpy.types.Object, mat=None):
        super(Cell, self).__init__(obj)

        # Set up effector collections
        self._effectors = bpy.data.collections.new(f"{obj.name}_effectors")
        bpy.context.scene.collection.children.link(self._effectors)
        # self.add_effector(ForceCollection.global_forces())  # link global forces

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
        """Name of the cell. Also defines the name of related forces and
        collections of effectors.
        """
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
        """Copies the cell.

        The underlying Blender object, object data, and material if applicable.

        Warning:
            Any settings that use custom collections will not be updated. It is
            advised that physics modifiers are set up again for the daughter
            cell after calling `copy()`.

        Returns:
            A new cell with copied object and mesh data.
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
    def obj_eval(self) -> bpy.types.ID.evaluated_get:
        """The evaluated object.

        Note:
            `Blender API Documentation > evaluated_get(depsgraph)`__

        __ https://docs.blender.org/api/current/bpy.types.ID.html?highlight=evaluated_get#bpy.types.ID.evaluated_get
        """
        dg = bpy.context.evaluated_depsgraph_get()
        obj_eval = self.obj.evaluated_get(dg)
        return obj_eval

    def vertices(self, local_coords: bool = False) -> list[Vector]:
        """Returns the vertices of the mesh representation of the cell.

        Args:
            local_coords: if `True`, coordinates are returned in local object space rather than world space.

        Returns:
            List of coordinates of vertices.
        """
        verts = self.obj_eval.data.vertices
        if local_coords:
            return [v.co for v in verts]
        else:
            return [self.obj_eval.matrix_world @ v.co for v in verts]

    def volume(self) -> float:
        """Calculates the volume of the cell.

        Returns:
            The signed volume of the cell (with physics evaluated).
        """
        bm = bmesh.new()
        bm.from_mesh(self.obj_eval.to_mesh())
        bm.transform(self.obj_eval.matrix_world)
        volume = bm.calc_volume()
        bm.free()

        return volume

    def COM(self, local_coords: bool = False) -> Vector:
        """Calculates the center of mass of a cell.

        Args:
            local_coords: if `True`, coordinates are returned in local object
                space rather than world space.

        Returns:
            The vector representing the center of mass of the cell.
        """
        vert_coords = self.vertices(local_coords)
        com = Vector(np.mean(vert_coords, axis=0))
        return com

    def _get_eigenvector(self, n: int) -> Axis:
        """Returns the nth eigenvector (axis) in object space as a line defined
        by two Vectors.

        This function calculates the nth eigenvector of the covariance matrix
        defined by the cell's vertices. The eigenvector axis is defined by two
        points: the vertices in the direction of the eigenvector with the
        smallest and largest projections.

        Args:
            n: The index of the eigenvector to return.

        Returns:
            An axis defined by the eigenvector and the vertices at the
            extreme projections along this vector.
        """
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
        """Returns the major axis of the cell."""
        return self._get_eigenvector(0)

    def minor_axis(self) -> Axis:
        """Returns the minor axis of the cell."""
        return self._get_eigenvector(1)

    def divide(self, division_logic) -> tuple["Cell", "Cell"]:
        """Cause the cell to divide into two daughter cells.

        This function causes the cell to divide into two daughter cells according
        to the provided division logic. The daughter cells inherit the cell type
        of the mother cell.

        Args:
            division_logic: The division logic to use, which handles the
                creation of two cells from the original cell.

        Returns:
            A tuple of two daughter cells, resulting from the division of the
            mother cell.
        """
        # TODO: rewrite code to make it clearer that there are two daughter cells splitting
        # from a mother cell.
        mother, daughter = division_logic.make_divide(self)
        if mother.celltype:
            mother.celltype.add_cell(daughter)
        return mother, daughter

    def recenter(self):
        """Recenter the cell origin to the center of mass of the cell."""
        bm = bmesh.new()
        bm.from_mesh(self.obj_eval.to_mesh())

        com = self.COM()
        bmesh.ops.translate(bm, verts=bm.verts, vec=-self.COM(local_coords=True))

        bm.to_mesh(self.obj.data)
        bm.free()

        self.loc = com

    def remesh(self, voxel_size: float = 0.25, smooth: bool = True):
        """Remesh the underlying mesh representation of the cell.

        Remeshing is done using the built-in `voxel_remesh()`.

        Args:
            voxel_size: The resolution used for the remesher (smaller means more
            polygons).  smooth: If true, the final cell faces will appear
            smooth.
        """
        # use of object ops is 2x faster than remeshing with modifiers
        self.obj.data.remesh_mode = "VOXEL"
        self.obj.data.remesh_voxel_size = voxel_size
        with bpy.context.temp_override(active_object=self.obj, object=self.obj):
            bpy.ops.object.voxel_remesh()

        for f in self.obj.data.polygons:
            f.use_smooth = smooth

    def recolor(self, color: tuple[float, float, float]):
        """Recolors the material of the cell.

        This function changes the diffuse color of the cell's material to the
        specified color while preserving the alpha value. If the material uses
        nodes, it also updates the 'Base Color' input of any nodes that have it.

        Args:
            color: A tuple (r, g, b) representing the new color to apply.
        """
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
        """Retrieves the first modifier of the specified type from the
        underlying object representation of the cell.

        Args:
            type: The type of the modifier to search for.

        Returns:
            The first modifier of the specified type if found, otherwise None.
        """
        return next((m for m in self.obj.modifiers if m.type == type), None)

    @property
    def cloth_mod(self) -> Optional[ClothModifier]:
        """The cloth modifier of the cell if it exists, otherwise None."""
        return self.get_modifier("CLOTH")

    @property
    def collision_mod(self) -> Optional[CollisionModifier]:
        """The collision modifier of the cell if it exists, otherwise None."""
        return self.get_modifier("COLLISION")

    @property
    def physics_enabled(self) -> bool:
        """Whether physics is enabled for this cell."""
        return self._physics_enabled

    def _update_cloth(self):
        """Update the cloth modifier is correctly set to be affected by forces
        acting upon the cell.
        """
        if self.cloth_mod:
            self.cloth_mod.settings.effector_weights.collection = self._effectors

    def setup_physics(self, physics_constructor: PhysicsConstructor):
        """Set up the physics properties for the cell.

        This function initializes the physics properties of the cell using the given
        physics constructor. It then marks the physics as enabled.

        Args:
            physics_constructor: A function or callable that sets up the physics
                properties for the object.
        """
        physics_constructor(self.obj)
        self._update_cloth()
        self._physics_enabled = True

    def enable_physics(self):
        """Enable the physics simulation for the cell.

        This function re-enables the physics simulation for the cell by recreating
        the modifier stack from stored settings, updating the cloth modifier,
        and enabling any adhesion forces.

        Raises:
            RuntimeError: If physics is already enabled.
        """
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
        """
        Disable the physics simulation for the cell.

        This function disables the physics simulation for the cell by storing the
        current modifier settings, removing all modifiers, and disabling any adhesion
        forces.

        Raises:
            RuntimeError: If physics is not enabled.
        """
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
        """Stiffness of the membrane of the cell."""
        return self.cloth_mod.settings.tension_stiffness

    @stiffness.setter
    def stiffness(self, stiffness: float):
        self.cloth_mod.settings.tension_stiffness = stiffness
        self.cloth_mod.settings.compression_stiffness = stiffness
        self.cloth_mod.settings.shear_stiffness = stiffness

    @property
    def pressure(self) -> float:
        """Internal pressure of the cell."""
        return self.cloth_mod.settings.uniform_pressure_force

    @pressure.setter
    def pressure(self, pressure: float):
        self.cloth_mod.settings.uniform_pressure_force = pressure

    # ----- FORCES -----
    def add_effector(self, force: Force | ForceCollection):
        """Add a force or a collection of forces that affects this cell.

        Args:
            force: The force or collection of forces to add.
        """
        if isinstance(force, Force):
            self._effectors.objects.link(force.obj)
        elif isinstance(force, ForceCollection):
            self._effectors.children.link(force.collection)

    def remove_effector(self, force: Force | ForceCollection):
        """Remove a force or a collection of forces that affects this cell.

        Args:
            force: The force or collection of forces to remove.
        """
        if isinstance(force, Force):
            self._effectors.objects.unlink(force.obj)
        elif isinstance(force, ForceCollection):
            self._effectors.children.unlink(force.collection)

    def link_adhesion_force(self, force: AdhesionForce):
        """Set this cell as the origin of an adhesion force.

        Args:
            force: The adhesion force which originates from this cell.
        """
        self._hetero_adhesions.append(force)

    @property
    def homo_adhesion(self) -> AdhesionForce:
        """Homotypic adhesion force of the cell."""
        return self._homo_adhesion

    @homo_adhesion.setter
    def homo_adhesion(self, force: AdhesionForce):
        self._homo_adhesion = force

    @property
    def adhesion_forces(self) -> list[AdhesionForce]:
        """Heterotypic adhesion forces of the cell."""
        return [self._homo_adhesion] + self._hetero_adhesions

    @property
    def motion_force(self) -> Force:
        """Motion force of the cell."""
        return self._motion_force

    @motion_force.setter
    def motion_force(self, force: MotionForce):
        if self.motion_force:
            self.remove_effector(self.motion_force)
        self.add_effector(force)
        self._motion_force = force

    def move_towards(self, dir: Vector):
        """Sets the motion force to move the cell in a specified direction.

        This function sets the motion force to move the cell in the specified
        direction. If the cell does not have an associated motion force, it raises
        a RuntimeError.

        Args:
            dir (Vector): The direction in which to set the motion force.

        Raises:
            RuntimeError: If the cell does not have an associated motion force.
        """
        if not self._motion_force:
            raise RuntimeError(
                f"Cell {self.name} does not have an associated motion force!"
            )
        motion_loc = self.loc + dir.normalized() * (2 + self.major_axis().length())
        self._motion_force.set_loc(motion_loc, self.loc)


def create_cell(
    name: str,
    loc: tuple,
    color: tuple = None,
    physics_constructor: PhysicsConstructor = None,
    physics_on: bool = True,
    **kwargs,
) -> Cell:
    """Creates a new cell using the default cell type.

    Args:
        name: The name of the cell.
        loc: The location of the cell.
        color: The color of the cell.
        physics_constructor: A generator of physics properties of the cell.
        physics_on: Whether to enable physics for the cell.
        **kwargs: keyword arguments passed to the Blender mesh generator function.

    Returns:
        The newly created cell.
    """
    cell = CellType.default_celltype().create_cell(
        name, loc, color, physics_constructor, physics_on, **kwargs
    )
    return cell


def store_settings(mod: bpy.types.bpy_struct) -> dict:
    """Store the settings of a Blender modifier in a dictionary.

    Args:
        mod: The Blender modifier.

    Returns:
        A dictionary with the stored settings.
    """

    settings = {}
    for p in mod.bl_rna.properties:
        id = p.identifier
        if not p.is_readonly:
            settings[id] = getattr(mod, id)
        elif id in ["settings", "collision_settings", "effector_weights"]:
            settings[id] = store_settings(getattr(mod, id))
    return settings


def declare_settings(mod: bpy.types.bpy_struct, settings: dict):
    """Recursively apply stored settings to a Blender modifier.

    Args:
        mod: The Blender modifier to which the settings are applied.
        settings: A dictionary containing the settings.
    """
    for id, setting in settings.items():
        if isinstance(setting, dict):
            declare_settings(getattr(mod, id), settings[id])
        else:
            setattr(mod, id, setting)


class CellType:
    """A cell type.

    Cell types are represented in Blender by collection. Cells of the same cell
    type interact through homotypic adhesion forces, and interact with cells of
    different cell types through heterotypic adhesion forces. Cell types have
    default mesh creation settings (including color) and default physics
    settings (including cloth, collision, homotypic adhesion forces, and motion
    forces).

    Args:
        name: The name of the cell type.
        physics_enabled: Whether cells of this cell type are responsive to physics.

    Attributes:
        homo_adhesion_strength (int): Default homotypic adhesion strength
            between cells of this type.
        motion_strength (int): Default motion strength of
            cells of this type.
    """

    _mesh_kwargs = {}
    _physics_constructor = PhysicsConstructor(
        SubsurfConstructor,
        ClothConstructor,
        CollisionConstructor,
        # RemeshConstructor,
    )
    color = (0.007, 0.021, 0.3)
    _default_celltype = None

    def __init__(self, name: str, physics_enabled: bool = True):
        self._homo_adhesions = ForceCollection(name)
        self._cells = set()

        self._physics_enabled = physics_enabled
        self._hetero_adhesions = {}

        self.homo_adhesion_strength: int = 2000
        self.motion_strength: int = 0

    @staticmethod
    def default_celltype() -> "CellType":
        """Get the default cell type."""
        if CellType._default_celltype is None:
            CellType._default_celltype = CellType("default")
        return CellType._default_celltype

    @property
    def name(self) -> str:
        """Name of the cell type."""
        return self._homo_adhesions.name

    def add_cell(self, cell: Cell):
        """Add a cell to the cell type, activating its physics and constructing
        appropriate forces.

        This method sets the cell type of the added cell to this cell type and
        activates physics for the cell. It constructs and adds homotypic
        adhesion forces, heterotypic adhesion forces, and motion forces defined
        for this cell type.

        Args:
            cell: The cell to add to the cell type.
        """
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
        motion.hide()

    # TODO: implement correctly
    def _remove_cell(self, cell: Cell):
        pass
        # cell.celltype = None
        # self._cells.remove(cell)

    def create_cell(
        self,
        name: str,
        loc: tuple,
        color: tuple = None,
        physics_constructor: PhysicsConstructor = None,
        physics_on: bool = True,
        **mesh_kwargs,
    ) -> Cell:
        """Creates a new cell, then adds it to this cell type.

        Args:
            name: The name of the cell.
            loc: The location of the cell.
            color: The color of the cell.
            physics_constructor: A generator of physics properties of the cell.
            physics_on: Whether to enable physics for the cell.
            **mesh_kwargs: keyword arguments passed to the Blender mesh
                generator function.

        Returns:
            The newly created cell.

        Note:
            `create_cell` does not directly activate physics, which is handled
            by the `add_cell` function. This is significant for division, which
            will use the default `CellType` settings to set up cell physics,
            rather than any custom settings of the dividing cell.
        """
        mesh_kwargs = dict(self.__class__._mesh_kwargs, **mesh_kwargs)
        if color is None:
            color = self.__class__.color

        obj = create_mesh(name, loc, mesh="icosphere", **mesh_kwargs)
        bpy.context.scene.collection.objects.link(obj)

        mat = create_material(f"{name}_material", color=color) if color else None
        cell = Cell(obj, mat)
        # cell.remesh()

        # enable physics for cell
        if physics_constructor is None:
            physics_constructor = self.__class__._physics_constructor
        if self._physics_enabled and physics_on:
            cell.setup_physics(physics_constructor)

        self.add_cell(cell)
        return cell

    # TODO: create cells in shape or in specified locations of specified cell types
    def create_cells(celltype, locs=None, shape=None, size=None):
        pass

    @property
    def cells(self) -> list[Cell]:
        """The list of cells associated with this cell type."""
        return list(self._cells)

    def set_hetero_adhesion(self, other_celltype: "CellType", strength: float):
        """Set the default strength of heterotypic adhesion forces between this
        cell type and another.

        Args:
            other_celltype: The other cell type in relation to which heterotypic
                adhesion forces will be defined.
            strength: The force of the adhesion forces.
        """
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
    """A cell type that is reduced in complexity for faster rendering."""

    color = None
    physics_constructor = PhysicsConstructor(
        ClothConstructor,
        CollisionConstructor,
    )


class YolkType(CellType):
    """A larger cell type used for modeling embryonic yolks."""

    mesh_kwargs = {
        "size": 10,
        "subdivisions": 4,
    }
    physics_constructor = PhysicsConstructor(
        YolkClothConstructor,
        CollisionConstructor,
    )
    color = (0.64, 0.64, 0.64)
