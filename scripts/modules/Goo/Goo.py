# Goo.py - This is the Goo library. It contains all helper functions for Goo
# Goo is licensed under BSDv2

import bpy, bmesh, mathutils, numpy as np

def calculate_volume(obj): # obj = bpy.context.object
    dg = bpy.context.evaluated_depsgraph_get()
    obj_eval = obj.evaluated_get(dg)
    mesh_from_eval = obj_eval.to_mesh()
    bm = bmesh.new()
    bm.from_mesh(mesh_from_eval)
    volume = bm.calc_volume(signed = True)
    return volume

def get_major_axis(obj): # obj = bpy.context.object
    dg = bpy.context.evaluated_depsgraph_get()
    obj_eval = obj.evaluated_get(dg)
    vertices = obj_eval.data.vertices
    vert_coords = [(obj_eval.matrix_world @ v.co) for v in vertices]
    vert_coords = np.asarray(vert_coords)

    x = vert_coords[:, 0]
    x = x - np.mean(x)
    y = vert_coords[:, 1]
    y = y - np.mean(y)
    z = vert_coords[:, 2]
    z = z - np.mean(z)

    new_coords = np.vstack([x, y, z])
    cov_matrix = np.cov(new_coords)
    eigenvals, eigenvecs = np.linalg.eig(cov_matrix)

    sort_indices = np.argsort(eigenvals)
    major_x, major_y, major_z = eigenvecs[:, sort_indices[-1]]
    major_axis = (major_x, major_y, major_z)
    return major_axis

def calculate_contact_area(obj): # obj = bpy.context.object
    obj = bpy.context.active_object
    dg = bpy.context.evaluated_depsgraph_get()
    obj_eval = obj.evaluated_get(dg)
    mesh_from_eval = obj_eval.to_mesh()
    bm = bmesh.new()
    bm.from_mesh(mesh_from_eval)

    contact_area = 0

    for edge in bm.edges:
        angle = edge.calc_face_angle()
        if angle < 0.01:
            for face in edge.link_faces:
                face.select = True
                
    for face in bm.faces:
        if face.select == True:
            contact_area += face.calc_area()

    return contact_area

# Repair hole: after cell division, add new face and re-triangulate mesh
# may want to consider
def repair_hole(obj): 
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.mode_set(mode = 'EDIT')
    bpy.ops.mesh.select_mode(type="EDGE") 
    bpy.ops.mesh.select_all(action = 'SELECT')
    bpy.ops.mesh.edge_face_add()
    bpy.ops.object.mode_set(mode = 'OBJECT')
    #bpy.ops.mesh.select_all(action = 'DESELECT')
    #bpy.ops.object.mode_set(mode = 'OBJECT')
    #bpy.ops.object.modifier_add(type='REMESH')
    #bpy.context.object.modifiers["Remesh"].voxel_size = 0.3
    #bpy.ops.object.modifier_apply(modifier="Remesh")
    #bpy.ops.object.modifier_add(type='TRIANGULATE')
    #bpy.ops.object.modifier_apply(modifier="Triangulate")

def divide(obj): # obj = bpy.context.object
    # Get Center of Mass
    mother_name = obj.name
    daughter_name = mother_name + ".001"
    bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_VOLUME', center='MEDIAN')
    COM = obj.location
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    major_axis = get_major_axis(obj)
    bpy.ops.mesh.bisect(plane_co = COM, plane_no = major_axis, use_fill=False, flip=False)
    bpy.ops.object.mode_set(mode = 'OBJECT')
    # Separate the object
    obj = bpy.data.objects[mother_name]
    new_vertices = obj.data.vertices
    new_vert_coords = [(obj.matrix_world @ v.co) for v in new_vertices]
    new_vert_coords = np.asarray(new_vert_coords)

    separation_indices = []
    for i in range(len(new_vert_coords)):
        distance = mathutils.geometry.distance_point_to_plane(new_vert_coords[i], COM, major_axis)
        if distance > -0.05:
            separation_indices.append(i)
    # obj = bpy.context.active_object
    bpy.ops.object.mode_set(mode = 'EDIT') 
    bpy.ops.mesh.select_mode(type="VERT")
    bpy.ops.mesh.select_all(action = 'DESELECT')
    bpy.ops.object.mode_set(mode = 'OBJECT')
    for index in separation_indices:
        obj.data.vertices[index].select = True
    bpy.ops.object.mode_set(mode = 'EDIT') 
    bpy.ops.mesh.separate(type='SELECTED')
    bpy.ops.object.mode_set(mode = 'OBJECT')
    bpy.data.objects[mother_name].select_set(False)
    bpy.data.objects[daughter_name].select_set(False)
    bpy.context.view_layer.objects.active = obj
    repair_hole(obj)
    bpy.data.objects[mother_name].select_set(False)
    bpy.context.object.name = mother_name + "0"
    bpy.context.view_layer.objects.active = bpy.data.objects[daughter_name]
    repair_hole(bpy.data.objects[daughter_name])
    bpy.data.objects[daughter_name].select_set(False)
    bpy.context.object.name = mother_name + "1"

def mitosis_handler(scene):
    num_cells = len(bpy.data.collections["Cells"].objects)
    for i in range(num_cells):
        bpy.context.view_layer.objects.active = bpy.data.collections["Cells"].objects[i]
        cell = bpy.context.view_layer.objects.active
        volume = calculate_volume(cell)
        if volume > .3:
            divide(cell)
  
def make_cell(cell):
    #mesh = bpy.data.meshes.new()
    #cell = bmesh.ops.create_icosphere(mesh, 2, 2.0, insert_matrix_here, calc_uv = True)
    #bpy.ops.mesh.primitive_ico_sphere_add(radius = cell.radius, enter_editmode = cell.enter_editmode, align = cell.align, location = cell.location, scale = cell.scale)
    bpy.ops.mesh.primitive_round_cube_add(change=True, radius=cell.radius, size= cell.size, arc_div= cell.arcdiv, lin_div=0, div_type='CORNERS', odd_axis_align=False, no_limit=False)
    bpy.ops.object.modifier_add(type = 'SUBSURF')
    bpy.context.object.modifiers["Subdivision"].levels = cell.subdiv
    #bpy.ops.object.modifier_add(type = 'SOFT_BODY')
    bpy.ops.object.modifier_add(type = 'CLOTH')
    bpy.context.object.modifiers["Cloth"].settings.bending_model = 'LINEAR'
    bpy.context.object.modifiers["Cloth"].settings.quality = 5
    bpy.context.object.modifiers["Cloth"].settings.time_scale = 1
    bpy.context.object.modifiers["Cloth"].settings.mass = cell.vertex_mass
    bpy.context.object.modifiers["Cloth"].settings.air_damping = cell.air_damping
    bpy.context.object.modifiers["Cloth"].settings.tension_stiffness = 15
    bpy.context.object.modifiers["Cloth"].settings.compression_stiffness = 15
    bpy.context.object.modifiers["Cloth"].settings.shear_stiffness = 5
    bpy.context.object.modifiers["Cloth"].settings.bending_stiffness = 0.5
    bpy.context.object.modifiers["Cloth"].settings.tension_damping = 5
    bpy.context.object.modifiers["Cloth"].settings.compression_damping = 5
    bpy.context.object.modifiers["Cloth"].settings.shear_damping = 5
    bpy.context.object.modifiers["Cloth"].settings.bending_damping = 0.5
    bpy.context.object.modifiers["Cloth"].settings.use_pressure = True
    bpy.context.object.modifiers["Cloth"].settings.uniform_pressure_force = 2.8
    bpy.context.object.modifiers["Cloth"].settings.pressure_factor = 1
    bpy.context.object.modifiers["Cloth"].settings.fluid_density = 0
    bpy.context.object.modifiers["Cloth"].collision_settings.use_collision = True
    bpy.context.object.modifiers["Cloth"].collision_settings.distance_min = 0.015
    bpy.context.object.modifiers["Cloth"].collision_settings.impulse_clamp = 0 
    """
    bpy.context.object.modifiers["Softbody"].settings.use_goal = cell.softbody_goal 
    #bpy.context.object.modifiers["Softbody"].settings.friction = cell.softbody_friction
    bpy.context.object.modifiers["Softbody"].settings.mass = cell.mass
    bpy.context.object.modifiers["Softbody"].settings.pull = cell.edges_pull
    bpy.context.object.modifiers["Softbody"].settings.bend = cell.edges_bend
    bpy.context.object.modifiers["Softbody"].settings.use_self_collision = cell.self_collision
    bpy.context.object.modifiers["Softbody"].settings.ball_stiff = cell.self_collision_stiffness
    bpy.context.object.modifiers["Softbody"].settings.ball_size = cell.ball_size
    bpy.context.object.modifiers["Softbody"].settings.friction = 0.2
    """
    bpy.ops.object.modifier_add(type='COLLISION')
    bpy.context.object.collision.use_culling = False
    #bpy.context.object.collision.damping = 0.579821
    #bpy.context.object.collision.thickness_outer = 0.02
   # bpy.context.object.collision.thickness_inner = 0.2
    #bpy.context.object.collision.cloth_friction = 5
    bpy.ops.object.forcefield_toggle()
    bpy.context.object.field.type = 'FORCE'
    #bpy.context.object.field.strength = -800
    bpy.context.object.field.strength = -1
    bpy.context.object.field.shape = 'SURFACE'
    bpy.context.object.name = cell.name
    
class Cell():
    def __init__(self, name_string, loc):
        self.name = name_string
        self.radius = 1
        self.enter_editmode = False
        self.align = 'WORLD'
        self.location = loc
        self.size = (2, 2, 2)
        self.scale = (1, 1, 1)
        self.arcdiv = 8
        self.subdiv = 2
        self.vertex_mass = 0.3
        self.density = 1.0
        #self.init_volume = ((4/3)*np.pi*(self.radius)**3)*(10**-12) # fix to be in cm^3
        self.init_volume = ((4/3)*np.pi*(self.radius)**3)
        #self.mass = (self.density*self.init_volume)**1000
        self.mass = self.density*self.init_volume
        self.edges_pull = 0.0
        self.edges_bend = 0.2
        self.air_damping = 10
        self.self_collision = True
        self.self_collision_stiffness = 0.01
        self.ball_size = 10
        self.softbody_goal = False
        self.sotbody_friction = 0.2
        
    def get_blender_object(self):
        obj = bpy.data.objects[self.name]
        return obj