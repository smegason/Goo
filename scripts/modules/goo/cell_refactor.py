import bpy
import numpy as np

from .utils import *
from .force import MotionForce, ForceCollection, create_adhesion, create_motion
from .growth import *
from .circuits import GeneRegulatoryNetwork as GRN, Circuit


# TODO: perhaps refactor forces into CellBody subobject
class Cell(BlenderObject):
    def __init__(self):
        self.name = ""
        self.obj = None
        self.mat = None
        self.celltype = None

        self.adhesion_forces = None
        self.motion_force: MotionForce = None
        self.effectors = []

        self.growth_controller: GrowthController = None
        self.grn: GRN = None

        self.just_divided = False
        self.physics_enabled = False

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    def copy(self):
        pass

    # ===== Physical calculations =====
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

        projections = verts @ axis
        min_index = np.argmin(projections)
        max_index = np.argmax(projections)
        first_vertex = Vector(verts[min_index])
        last_vertex = Vector(verts[max_index])

        return Axis(Vector(axis), first_vertex, last_vertex, self.obj_eval.matrix_world)

    def major_axis(self) -> Axis:
        """Returns the major axis of the cell."""
        return self._get_eigenvector(0)

    def minor_axis(self) -> Axis:
        """Returns the minor axis of the cell."""
        return self._get_eigenvector(1)

    def radius(self):
        """Calculate the radius based on the major and minor axes of the cell."""
        major_axis = self.major_axis()
        minor_axis = self.minor_axis()
        return (major_axis.length() + minor_axis.length()) / 2

    # ===== Mesh operations =====
    @property
    def obj_eval(self) -> bpy.types.ID:
        """The evaluated object.

        Note:
            See the `Blender API Documentation for evaluated_get(depsgraph)
            <https://docs.blender.org/api/current/bpy.types.ID.html?highlight=evaluated_get#bpy.types.ID.evaluated_get>`__.
        """
        dg = bpy.context.evaluated_depsgraph_get()
        obj_eval = self.obj.evaluated_get(dg)
        return obj_eval

    def vertices(self, local_coords: bool = False) -> list[Vector]:
        """Returns the vertices of the mesh representation of the cell.

        Args:
            local_coords: if `True`, coordinates are returned in local object space
            rather than world space.

        Returns:
            List of coordinates of vertices.
        """
        verts = self.obj_eval.data.vertices
        if local_coords:
            return [v.co.copy() for v in verts]
        else:
            matrix_world = self.obj_eval.matrix_world
            return [matrix_world @ v.co for v in verts]

    def recenter(self, origin=True, forces=True):
        """Recenter cell origin to center of mass of cell, and center forces to that same origin."""
        pass

    def remesh(self, voxel_size: float = 0.55, smooth: bool = False):
        pass

    def recolor(self, color: tuple[float, float, float]):
        """Recolors the material of the cell.

        This function changes the diffuse color of the cell's material to the
        specified color while preserving the alpha value. If the material uses
        nodes, it also updates the 'Base Color' input of any nodes that have it.

        Args:
            color: A tuple (r, g, b) representing the new color to apply.
        """
        r, g, b = color
        _, _, _, a = self.mat.diffuse_color
        self.mat.diffuse_color = (r, g, b, a)

        if self.mat.use_nodes:
            for node in self.mat.node_tree.nodes:
                if "Base Color" in node.inputs:
                    _, _, _, a = node.inputs["Base Color"].default_value
                    node.inputs["Base Color"].default_value = r, g, b, a

    # ===== Physics operations =====
    def enable_physics(self):
        pass

    def disable_physics(self):
        pass

    @property
    def cloth_mod(self):
        pass

    @property
    def stiffness(self):
        pass

    @stiffness.setter
    def stiffness(self, stiffness: float):
        pass

    @property
    def pressure(self):
        pass

    @pressure.setter
    def pressure(self, pressure: float):
        pass

    def add_effector(self, force):
        pass

    # TODO: overload add_force to deal with homo, hetero, motion forces.
    def add_force(self, force):
        pass

    # ===== Cell actions =====
    # TODO: integrate growth controller with division
    def divide(self, division_handler):
        self.just_divided = True
        self.growth_controller.step_divided(self.volume())
        pass

    def move_towards(self, direction: Vector, strength=None):
        motion_loc = self.loc + dir.normalized() * (2 + self.radius())
        self.motion_force.set_loc(motion_loc, self.loc)
        if strength is not None:
            self.motion_force.strength = strength

    def step_growth(self, dt=1):
        if self.just_divided:
            self.just_divided = False
        else:
            new_pressure = self.growth_controller.step_growth(self.volume(), dt)
            self.pressure = new_pressure

    def step_grn(self, dt=1):
        """Calculate and update the metabolite concentrations of the gene regulatory network after 1 time step."""
        com = self.COM()
        radius = self.radius()

        self.grn.update_concs(center=com, radius=radius, dt=dt)
        self["genes"] = {str(k): v for k, v in self.metabolites.items()}

    # ===== Gene regulatory network items =====
    @property
    def metabolites(self):
        return self.grn.concs

    @metabolites.setter
    def metabolites(self, metabolites):
        self["genes"] = {str(k): v for k, v in metabolites.items()}
        self.grn.concs = metabolites


class DiffusionSystem:
    def __init__(self):
        pass

    def get_concentration(self, voxel):
        pass

    def get_total_concentration(self, center, radius):
        pass


class CellType:
    """A director class that takes a cell pattern as an instruction set to create
    new cells."""

    _default_celltype = None

    def __init__(
        self,
        name,
        pattern="simple",
        homo_adhesion_strength=2000,
        hetero_adhesion_strengths={},
        motion_strength=0,
    ):
        self.name = name

        if type(pattern) == str:
            match pattern:
                case "simple":
                    self.pattern = SimplePattern()
                case "yolk":
                    self.pattern = YolkPattern()
                case _:
                    self.pattern = StandardPattern()
        else:
            self.pattern = pattern

        self.homo_adhesion_strength = homo_adhesion_strength
        self._homo_adhesion_collection = ForceCollection(name)

        self.motion_strength = motion_strength

        self.hetero_adhesion_strengths = {}
        self._hetero_adhesion_collections = {}
        for celltype, strength in hetero_adhesion_strengths.items():
            self.set_hetero_adhesion_strength(celltype, strength)

    @staticmethod
    def default_celltype() -> "CellType":
        """Get the default cell type."""
        if CellType._default_celltype is None:
            CellType._default_celltype = CellType("default")
        return CellType._default_celltype

    def create_cell(
        self,
        name,
        loc,
        color: tuple = None,
        mesh_kwargs: dict = {},
        physics_enabled: bool = True,
        physics_constructor: PhysicsConstructor = None,
        growth_enabled: bool = True,
        growth_type: Growth = Growth.LINEAR,
        growth_rate: float = 1,
        target_volume: float = 30,
        genes_enabled: bool = True,
        circuits: list[Circuit] = None,
        metabolites: dict[str, float] = {},
    ):
        self.pattern.build_obj(
            name,
            loc,
            color=color,
            **mesh_kwargs,
        )
        if physics_enabled:
            self.pattern.build_physics(
                physics_constructor=physics_constructor,
            )
            self.pattern.build_forces(
                self.homo_adhesion_strength,
                self._homo_adhesion_collection,
                self.motion_strength,
                self.hetero_adhesion_strengths,
                self._hetero_adhesion_collections,
            )
        if growth_enabled:
            self.pattern.build_growth_controller(
                growth_type=growth_type,
                growth_rate=growth_rate,
                target_volume=target_volume,
            )
        if genes_enabled:
            self.pattern.build_network(
                circuits=circuits,
                metabolites=metabolites,
            )
        cell = self.pattern.retrieve_cell()
        cell.celltype = self
        return cell

    def set_hetero_adhesion_strength(self, other_celltype: "CellType", strength: float):
        outgoing_collections = ForceCollection(f"{self.name}_to_{other_celltype.name}")
        incoming_collections = ForceCollection(f"{other_celltype.name}_to_{self.name}")

        self.hetero_adhesion_strengths[other_celltype] = strength
        self._hetero_adhesion_collections[other_celltype] = (
            outgoing_collections,
            incoming_collections,
        )

        self.hetero_adhesion_strengths[self] = strength
        other_celltype._hetero_adhesion_collections[self] = (
            incoming_collections,
            outgoing_collections,
        )


class CellPattern:
    """Builders of different cell parts. The CellType class takes these
    patterns in order to create new cells."""

    mesh_kwargs = {}
    color = None
    physics_constructor = None

    circuits = []
    metabolites = {}

    def __init__(self):
        self.reset()

    def reset(self):
        self._cell = Cell()

    def retrieve_cell(self):
        """Retrieves constructed cell, links cell to the scene, and resets the director."""
        # Link cell to scene
        cell = self._cell
        bpy.context.scene.collection.objects.link(cell.obj)

        self.reset()
        return cell

    @staticmethod
    def _override(base_option, user_option):
        if user_option is None:
            return base_option
        if isinstance(base_option, dict):
            return dict(base_option, **user_option)
        return user_option

    def build_obj(
        self,
        name,
        loc,
        color=None,
        **mesh_kwargs,
    ):
        self._cell.name = name
        mesh_kwargs = self._override(self.__class__.mesh_kwargs, mesh_kwargs)
        color = self._override(self.__class__.color, color)

        # Create cell object
        obj = create_mesh(name, loc, mesh="icosphere", **mesh_kwargs)
        if color is not None:
            mat = create_material(f"{name}_material", color=color)
            obj.data.materials.append(mat)
        self._cell.obj = obj
        self._cell.mat = mat

    def build_physics(
        self,
        physics_constructor=None,
    ):
        # Add physics modifiers to cell object
        physics_constructor = self._override(
            self.__class__.physics_constructor, physics_constructor
        )
        if physics_constructor is not None:
            physics_constructor(self._cell)
            self._cell.enable_physics()

    def build_forces(
        self,
        homo_adhesion_strength: int,
        homo_adhesion_collection: ForceCollection,
        motion_strength: int,
        hetero_adhesion_strengths: dict[CellType, int],
        hetero_adhesion_collections: dict[CellType, tuple[ForceCollection]],
    ):
        homo_adhesion = create_adhesion(homo_adhesion_strength, obj=self._cell.obj)
        homo_adhesion_collection.add_force(homo_adhesion)

        self._cell.add_force(homo_adhesion)
        self._cell.add_effector(homo_adhesion_collection)

        for celltype in hetero_adhesion_strengths.keys():
            strength = hetero_adhesion_strengths[celltype]
            incoming, outgoing = hetero_adhesion_collections[celltype]

            hetero_adhesion = create_adhesion(
                strength,
                name=f"{self._cell.name}_to_{celltype.name}",
                loc=self._cell.loc,
            )
            outgoing.add_force(hetero_adhesion)

            self._cell.add_force(hetero_adhesion)
            self._cell.add_effector(incoming)

        motion_force = create_motion(
            name=f"{self._cell.name}_motion",
            loc=self._cell.loc,
            strength=motion_strength,
        )
        self._cell.motion_force = motion_force
        motion_force.hide()

    def build_growth_controller(
        self, growth_type: Growth, growth_rate: float = 1, target_volume=30
    ):
        controller = PIDController(
            self._cell.volume(),
            growth_type=growth_type,
            growth_rate=growth_rate,
            target_volume=target_volume,
        )
        self._cell.growth_controller = controller

    def build_network(self, circuits=None, metabolites=None):
        circuits = self._override(self.__class__.circuits, circuits)
        metabolites = self._override(self.__class__.metabolites, metabolites)

        grn = GRN()
        grn.load_circuits(*circuits)
        self._cell.grn = grn
        self._cell.metabolites = metabolites


class StandardPattern(CellPattern):
    mesh_kwargs = {}
    physics_constructor = PhysicsConstructor(
        SubsurfConstructor,
        ClothConstructor,
        CollisionConstructor,
        RemeshConstructor,
    )
    color = (0.007, 0.021, 0.3)


class SimplePattern(CellPattern):
    """A cell type that is reduced in complexity for faster rendering."""

    color = None
    physics_constructor = PhysicsConstructor(
        ClothConstructor,
        CollisionConstructor,
    )


class YolkPattern(CellPattern):
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
