import bpy


class Cell:
    def __init__(self, name, loc):
        self.name = name
        self.obj = self._create_obj(loc)

    def set_celltype(self, celltype):
        self.celltype = celltype

    def link_circuit(self, circuit):
        self.circuit = circuit

    # --- Blender Functions ---
    def _create_obj(self, loc, radius=1, scale=(1, 1, 1), rotation=(0, 0, 0)):
        # creates blender object and selects it
        # returns: blender.type.Object object

        bpy.ops.mesh.primitive_ico_sphere_add(
            subdivisions=2,
            radius=radius,
            calc_uvs=True,
            align="WORLD",
            location=loc,
            rotation=rotation,
            scale=scale,
        )

        obj = bpy.context.object
        obj.name = self.name
        return obj

    @staticmethod
    def test():
        print("test")


class CellType:
    def __init__(self):
        pass

    def add_cell(self, cell):
        cell.set_celltype(self)
        self.cells.append(cell)


def create_cell(celltype, loc):
    cell = Cell(loc)
    celltype.add_cell(cell)


def create_cells(celltype, locs=None, shape=None, size=None):
    # create cells in shape or in specified locations of specified cell types
    pass
