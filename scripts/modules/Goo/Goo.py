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

#def get_division_angles(axis, prev_axis):
def get_division_angles(axis):
    z_axis = np.array((0, 0, 1))
    division_axis = axis/np.linalg.norm(axis)
    dot_product = np.dot(z_axis, division_axis)
    phi = np.arccos(dot_product)

    proj_axis_on_z = dot_product*z_axis
    proj_xy_axis = division_axis - proj_axis_on_z
    proj_xy_axis = proj_xy_axis/np.linalg.norm(proj_xy_axis)
    dot_product = np.dot((1, 0, 0), proj_xy_axis)
    theta = np.arccos(dot_product)

    #prev_axis_norm = prev_axis/np.linalg.norm(prev_axis)
    #dot_product = np.dot(division_axis, prev_axis_norm)
    #div_angle_change = np.arccos(dot_product)
    #return phi, theta, div_angle_change
    print("GOT DIVISION ANGLES")
    return phi, theta

def calculate_contact_area(obj): # obj = bpy.context.object
    dg = bpy.context.evaluated_depsgraph_get()
    obj_eval = obj.evaluated_get(dg)
    mesh_from_eval = obj_eval.to_mesh()
    bm = bmesh.new()
    bm.from_mesh(mesh_from_eval)

    contact_area = 0
    for face in bm.faces:
        face.select = False

    for edge in bm.edges:
        angle = edge.calc_face_angle_signed()
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
    #bpy.context.view_layer.objects.active = obj
    bpy.ops.object.mode_set(mode = 'EDIT')
    bpy.ops.mesh.select_mode(type="EDGE") 
    bpy.ops.mesh.select_all(action = 'SELECT')
    bpy.ops.mesh.edge_face_add()
    bpy.ops.mesh.select_all(action = 'DESELECT')
    bpy.ops.object.mode_set(mode = 'OBJECT')
    bpy.ops.object.modifier_add(type='REMESH')
    bpy.context.object.modifiers["Remesh"].voxel_size = 0.2
    bpy.ops.object.modifier_apply(modifier="Remesh")

#def divide(obj, tree): # obj = bpy.context.object
def sep(obj):
    # Get Center of Mass
    bpy.ops.object.mode_set(mode='OBJECT')
    dg = bpy.context.evaluated_depsgraph_get()
    obj = obj.evaluated_get(dg)
    mother_name = obj.name
    daughter_name = mother_name + ".001"
    bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS')
    COM = obj.location
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    # Get Major Axis
    major_axis = get_major_axis(obj)
    #phi, theta = get_division_angles(major_axis)
    #print(mother_name + " PHI: " + str(phi) + ", THETA: " + str(theta))
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
    bpy.data.objects[mother_name].select_set(True)
    return mother_name, daughter_name, COM, major_axis

def select_translate_cell(obj, COM, major_axis):
    com = np.copy(COM)
    print("OLD COM: " + str(COM))
    bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS')
    new_COM = obj.location
    print("NEW COM: " + str(new_COM))
    print(com)
    distance = mathutils.geometry.distance_point_to_plane(new_COM, com, major_axis)
    print(distance)
    if distance > 0:
        return 0
    else:
        return 1

def translate_cell(obj, major_axis):
    magnitude = ((major_axis[0])**2 + (major_axis[1])**2 + (major_axis[2])**2)**(0.5)
    new_vec = (0.1*major_axis[0]/magnitude, 0.1*major_axis[1]/magnitude, 0.1*major_axis[2]/magnitude)
    bpy.ops.transform.translate(value = new_vec)

def divide(obj):
    obj.select_set(True)
    m_name, d_name, COM, major_axis = sep(obj)
    print(major_axis)
    val = select_translate_cell(bpy.data.objects[m_name], COM, major_axis)
    if val == 0:
        translate_cell(bpy.data.objects[m_name], major_axis)
        bpy.data.objects[m_name].select_set(False)
    else:
        bpy.data.objects[d_name].select_set(True)
        bpy.data.objects[m_name].select_set(False)
        translate_cell(bpy.data.objects[d_name], major_axis)
        bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS')
        bpy.data.objects[d_name].select_set(False)
    bpy.data.objects[m_name].select_set(True)
    bpy.context.view_layer.objects.active = bpy.data.objects[m_name]
    repair_hole(bpy.data.objects[m_name])
    bpy.context.object.name = m_name + "0"
    daughter1 = Cell(m_name + "0", loc = bpy.context.object.location)
    daughter1.data['mother'] = m_name
    daughter1.data['daughters'] = ['none', 'none']
    bpy.data.objects[m_name + "0"].select_set(False)
    bpy.data.objects[d_name].select_set(True)
    bpy.context.view_layer.objects.active = bpy.data.objects[d_name]
    repair_hole(bpy.data.objects[d_name])
    bpy.context.object.name = m_name + "1"
    bpy.data.objects[m_name + "1"].select_set(False)
    daughter2 = Cell(bpy.context.object.name, loc = bpy.context.object.location)
    daughter2.data['mother'] = m_name
    daughter2.data['daughters'] = ['none', 'none']
    return daughter1, daughter2
    
def make_mesh(cell):
    try:
        bpy.ops.mesh.primitive_round_cube_add(change = False, radius=cell.data['radius'], size= cell.data['size'], arc_div= cell.data['arcdiv'], lin_div=0, div_type='CORNERS', odd_axis_align=False, no_limit=False, location = cell.data['location'])
    except:
        print("To enable RoundCube creation for Cells you must go to Edit->Preferences->AddOns->Add Mesh:ExtraObjects and check the box to enable it")
    #bpy.ops.mesh.primitive_cube_add(size= cell.size, location = cell.location, align = 'WORLD', scale = cell.scale)
    bpy.context.object.name = cell.data['name']
    bpy.context.view_layer.objects.active = bpy.data.objects[cell.data['name']]
    bpy.ops.object.modifier_add(type = 'SUBSURF')
    bpy.context.object.modifiers["Subdivision"].levels = cell.data['subdiv']
    
def turn_off_physics():
    bpy.ops.object.modifier_remove(modifier="Cloth")
    
def turn_on_physics():
    bpy.ops.object.modifier_add(type = 'CLOTH')
    bpy.context.object.modifiers["Cloth"].settings.bending_model = 'LINEAR'
    bpy.context.object.modifiers["Cloth"].settings.quality = 5
    bpy.context.object.modifiers["Cloth"].settings.time_scale = 1
    bpy.context.object.modifiers["Cloth"].settings.mass = 0.3
    bpy.context.object.modifiers["Cloth"].settings.air_damping = 10
    bpy.context.object.modifiers["Cloth"].settings.tension_stiffness = 15
    bpy.context.object.modifiers["Cloth"].settings.compression_stiffness = 15
    bpy.context.object.modifiers["Cloth"].settings.shear_stiffness = 5
    bpy.context.object.modifiers["Cloth"].settings.bending_stiffness = 15
    bpy.context.object.modifiers["Cloth"].settings.tension_damping = 5
    bpy.context.object.modifiers["Cloth"].settings.compression_damping = 5
    bpy.context.object.modifiers["Cloth"].settings.shear_damping = 5
    bpy.context.object.modifiers["Cloth"].settings.bending_damping = 0.5
    bpy.context.object.modifiers["Cloth"].settings.use_pressure = True
    bpy.context.object.modifiers["Cloth"].settings.uniform_pressure_force = 2.8
    bpy.context.object.modifiers["Cloth"].settings.pressure_factor = 1
    bpy.context.object.modifiers["Cloth"].settings.fluid_density = 1
    bpy.context.object.modifiers["Cloth"].collision_settings.use_collision = True
    bpy.context.object.modifiers["Cloth"].collision_settings.distance_min = 0.015
    bpy.context.object.modifiers["Cloth"].collision_settings.impulse_clamp = 0 
     

#def mitosis_handler(scene, tree):
def mitosis_handler(scene):
    num_cells = len(bpy.data.collections["Cells"].objects)
    for i in range(num_cells):
        bpy.context.view_layer.objects.active = bpy.data.collections["Cells"].objects[i]
        cell = bpy.context.view_layer.objects.active
        volume = calculate_volume(cell)
        if volume > .3:
            cell_tree = divide(cell, cell_tree)
    cell_tree.show()

def div_handler(scene):
    num_cells = len(bpy.data.collections["Cells"].objects)
    current_frame = bpy.data.scenes[0].frame_current
    if current_frame % 30 == 0:
        print("DIVIDING CELLS")
        for i in range(num_cells):
            cell_name = bpy.data.collections["Cells"].objects[i].name
            cell = bpy.data.objects[cell_name]
            bpy.data.objects[cell_name].select_set(True)
            bpy.context.view_layer.objects.active = bpy.data.objects[cell_name]
            print(cell_name)
            turn_off_physics()
            d1, d2 = divide(cell)
            bpy.data.objects[d1.data["name"]].select_set(True)
            bpy.data.objects[d2.data["name"]].select_set(False)
            bpy.context.view_layer.objects.active = bpy.data.objects[d1.data["name"]]
            turn_on_physics()
            bpy.data.objects[d1.data["name"]].select_set(False)
            bpy.data.objects[d2.data["name"]].select_set(True)
            bpy.context.view_layer.objects.active = bpy.data.objects[d2.data["name"]]
            turn_on_physics()
            bpy.data.objects[d2.data["name"]].select_set(False)

def make_cell(cell):
    print ("make cell");
    #mesh = bpy.data.meshes.new()
    #cell = bmesh.ops.create_icosphere(mesh, 2, 2.0, insert_matrix_here, calc_uv = True)
    #bpy.ops.mesh.primitive_ico_sphere_add(radius = cell.radius, enter_editmode = cell.enter_editmode, align = cell.align, location = cell.location, scale = cell.scale)
    bpy.ops.mesh.primitive_round_cube_add(change = False, radius=cell.data['radius'], size= cell.data['size'], arc_div= cell.data['arcdiv'], lin_div=0, div_type='CORNERS', odd_axis_align=False, no_limit=False, location = cell.data['location'])
    #bpy.ops.mesh.primitive_cube_add(size= cell.size, location = cell.location, align = 'WORLD', scale = cell.scale)
    bpy.context.object.name = cell.data['name']
    bpy.context.view_layer.objects.active = bpy.data.objects[cell.data['name']]
    bpy.ops.object.modifier_add(type = 'SUBSURF')
    bpy.context.object.modifiers["Subdivision"].levels = cell.data['subdiv']
    #bpy.ops.object.modifier_add(type = 'SOFT_BODY')
    bpy.ops.object.modifier_add(type = 'CLOTH')
    bpy.context.object.modifiers["Cloth"].settings.bending_model = 'LINEAR'
    bpy.context.object.modifiers["Cloth"].settings.quality = 5
    bpy.context.object.modifiers["Cloth"].settings.time_scale = 1
    bpy.context.object.modifiers["Cloth"].settings.mass = cell.data['vertex_mass']
    bpy.context.object.modifiers["Cloth"].settings.air_damping = cell.data['air_damping']
    bpy.context.object.modifiers["Cloth"].settings.tension_stiffness = 15
    bpy.context.object.modifiers["Cloth"].settings.compression_stiffness = 15
    bpy.context.object.modifiers["Cloth"].settings.shear_stiffness = 5
    bpy.context.object.modifiers["Cloth"].settings.bending_stiffness = 15
    bpy.context.object.modifiers["Cloth"].settings.tension_damping = 5
    bpy.context.object.modifiers["Cloth"].settings.compression_damping = 5
    bpy.context.object.modifiers["Cloth"].settings.shear_damping = 5
    bpy.context.object.modifiers["Cloth"].settings.bending_damping = 0.5
    bpy.context.object.modifiers["Cloth"].settings.use_pressure = True
    bpy.context.object.modifiers["Cloth"].settings.uniform_pressure_force = 2.8
    bpy.context.object.modifiers["Cloth"].settings.pressure_factor = 1
    bpy.context.object.modifiers["Cloth"].settings.fluid_density = 1
    bpy.context.object.modifiers["Cloth"].collision_settings.use_collision = True
    bpy.context.object.modifiers["Cloth"].collision_settings.distance_min = 0.015
    bpy.context.object.modifiers["Cloth"].collision_settings.impulse_clamp = 0 
    bpy.ops.object.modifier_add(type='COLLISION')
    bpy.context.object.collision.use_culling = False
    #bpy.context.object.collision.damping = 0.579821
    #bpy.context.object.collision.thickness_outer = 0.02
   # bpy.context.object.collision.thickness_inner = 0.2
    #bpy.context.object.collision.cloth_friction = 5
    #bpy.ops.object.forcefield_toggle()
    #bpy.context.object.field.type = 'FORCE'
    #bpy.context.object.field.strength = -600
    #bpy.context.object.field.strength = 0
    #bpy.context.object.field.shape = 'POINT'
    #bpy.context.object.name = cell.name
    if (bpy.data.materials.get(cell.data['material'])):
        bpy.context.active_object.data.materials.append(bpy.data.materials.get(cell.data['material']))
    else:
        print ("No material ", cell.data['material'])
        
class Cell():
    def __init__(self, name_string, loc, material=""):
        self.data = {
            'ID': 0,
            'name': name_string,
            'radius': 1,
            'enter_editmode': False,
            'align': 'WORLD',
            'location': loc,
            'material': material,
            'size': (1, 1, 1),
            'scale': (1, 1, 1),
            'arcdiv': 8,
            'subdiv': 2,
            'vertex_mass': 0.3,
            'density': 1.0,
            'edges_pull': 0.0,
            'edges_bend': 0.2,
            'air_damping': 10,
            'self_collision': True,
            'self_collision_stiffness': 0.01,
            'ball_size': 10,
            'softbody_goal': False,
            'sotbody_friction': 0.2,
            'phi': 0,
            'theta': 0,
            'div_axis': (0, 0, 0),
            'mother': 'none',
            'daughters': ['none', 'none']
        }
        self.data['init_volume'] = ((4/3)*np.pi*(self.data['radius'])**3)
        self.data['mass'] = self.data['density']*self.data['init_volume']
        """ self.id = 0
        self.name = name_string
        self.radius = 1
        self.enter_editmode = False
        self.align = 'WORLD'
        self.location = loc
        self.size = (1, 1, 1)
        #self.size = 2.0
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
        self.phi = 0
        self.theta = 0
        self.div_axis = (0, 0, 0) """
    def get_blender_object(self):
        obj = bpy.data.objects[self.data["name"]]
        return obj

    def divide(self):
        obj = self.get_blender_object()
        obj.select_set(True)
        m_name, d_name, COM, major_axis = sep(obj)
        print(major_axis)
        val = select_translate_cell(bpy.data.objects[m_name], COM, major_axis)
        if val == 0:
            translate_cell(bpy.data.objects[m_name], major_axis)
            bpy.data.objects[m_name].select_set(False)
        else:
            bpy.data.objects[d_name].select_set(True)
            bpy.data.objects[m_name].select_set(False)
            translate_cell(bpy.data.objects[d_name], major_axis)
            bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS')
            bpy.data.objects[d_name].select_set(False)
        bpy.data.objects[m_name].select_set(True)
        bpy.context.view_layer.objects.active = bpy.data.objects[m_name]
        repair_hole(bpy.data.objects[m_name])
        bpy.context.object.name = m_name + "0"
        daughter1 = Cell(m_name + "0", loc = bpy.context.object.location)
        daughter1.data['mother'] = m_name
        daughter1.data['daughters'] = ['none', 'none']
        bpy.data.objects[m_name + "0"].select_set(False)
        bpy.data.objects[d_name].select_set(True)
        bpy.context.view_layer.objects.active = bpy.data.objects[d_name]
        repair_hole(bpy.data.objects[d_name])
        bpy.context.object.name = m_name + "1"
        bpy.data.objects[m_name + "1"].select_set(False)
        daughter2 = Cell(bpy.context.object.name, loc = bpy.context.object.location)
        daughter2.data['mother'] = m_name
        daughter2.data['daughters'] = ['none', 'none']
        return daughter1, daughter2

def add_material(mat):
    bpy.data.materials.new(name = mat.name)
    bpy.data.materials[mat.name].node_tree.nodes["Mix Shader"].inputs[0].default_value = mat.BDSF_1_fac
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[0].default_value = mat.BDSF_1_color
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[1].default_value = mat.BDSF_1_subsurface
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[2].default_value[0] = mat.BDSF_1_subsurf_radius[0]
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[2].default_value[1] = mat.BDSF_1_subsurf_radius[1]
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[2].default_value[2] = mat.BDSF_1_subsurf_radius[2]
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[3].default_value = mat.BDSF_1_subsurf_color
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[4].default_value = mat.BDSF_1_metallic
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[5].default_value = mat.BDSF_1_specular
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[6].default_value = mat.BDSF_1_specular_tint
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[7].default_value = mat.BDSF_1_roughness
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[8].default_value = mat.BDSF_1_anisotropic
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[9].default_value = mat.BDSF_1_anisotropic_rot
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[10].default_value = mat.BDSF_1_sheen
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[11].default_value = mat.BDSF_1_sheen_tint
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[12].default_value = mat.BDSF_1_clearcoat
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[13].default_value = mat.BDSF_1_clear_rough
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[14].default_value = mat.BDSF_1_IOR
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[15].default_value = mat.BDSF_1_transmission
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[16].default_value = mat.BDSF_1_transmission_rough
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[17].default_value = mat.BDSF_1_emission_color
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[18].default_value = mat.BDSF_1_emission_strength
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[19].default_value = mat.BDSF_1_alpha
    bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[0].default_value = mat.BDSF_2_hue
    bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[1].default_value = mat.BDSF_2_saturation
    bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[2].default_value = mat.BDSF_2_value
    bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[3].default_value = mat.BDSF_2_fac
    bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[4].default_value = mat.BDSF_2_color
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[1].default_value = mat.BDSF_2_subsurface
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[2].default_value[0] = mat.BDSF_2_subsurf_radius[0]
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[2].default_value[1] = mat.BDSF_2_subsurf_radius[1]
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[2].default_value[2] = mat.BDSF_2_subsurf_radius[2]
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[3].default_value = mat.BDSF_2_subsurf_color
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[4].default_value = mat.BDSF_2_metallic
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[5].default_value = mat.BDSF_1_specular
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[6].default_value = mat.BDSF_2_specular_tint
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[7].default_value = mat.BDSF_2_roughness
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[8].default_value = mat.BDSF_2_anisotropic
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[9].default_value = mat.BDSF_2_anisotropic_rot
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[10].default_value = mat.BDSF_2_sheen
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[11].default_value = mat.BDSF_2_sheen_tint
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[12].default_value = mat.BDSF_2_clearcoat
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[13].default_value = mat.BDSF_2_clear_rough
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[14].default_value = mat.BDSF_2_IOR
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[15].default_value = mat.BDSF_2_transmission
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[16].default_value = mat.BDSF_2_transmission_rough
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[17].default_value = mat.BDSF_2_emission_color
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[18].default_value = mat.BDSF_2_emission_strength
    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[19].default_value = mat.BDSF_2_alpha
    bpy.context.object.active_material.blend_method = mat.blend_method
    bpy.context.object.active_material.shadow_method = mat.shadow_method
    bpy.context.object.active_material.refraction_depth = mat.refraction_depth
    bpy.context.object.active_material.pass_index = mat.pass_index
    bpy.context.object.active_material.diffuse_color = mat.diffuse_color
    bpy.context.object.active_material.metallic = mat.metallic
    bpy.context.object.active_material.roughness = mat.roughness
    bpy.ops.node.add_node(type="ShaderNodeTexNoise", use_transform=True)
    bpy.ops.node.select(wait_to_deselect_others=True, mouse_x=114, mouse_y=765, extend=False, deselect_all=True)
    bpy.ops.node.link(detach=False)

def add_material2(mat_name):
    mat = bpy.data.materials.get(mat_name)
    if mat is None:
        print ("Material not found. Make new")
        mat = bpy.data.materials.new(name=mat_name)
     
    if mat.node_tree:
        mat.node_tree.links.clear()
        mat.node_tree.nodes.clear()
            
    bpy.data.materials[mat.name].use_nodes = True
    
    #soap bubble shader
#    bpy.data.materials[mat.name].node_tree.nodes["Mix Shader"].inputs[0].default_value = mat.BDSF_1_fac
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[0].default_value = mat.BDSF_1_color
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[1].default_value = mat.BDSF_1_subsurface
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[2].default_value[0] = mat.BDSF_1_subsurf_radius[0]
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[2].default_value[1] = mat.BDSF_1_subsurf_radius[1]
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[2].default_value[2] = mat.BDSF_1_subsurf_radius[2]
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[3].default_value = mat.BDSF_1_subsurf_color
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[4].default_value = mat.BDSF_1_metallic
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[5].default_value = mat.BDSF_1_specular
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[6].default_value = mat.BDSF_1_specular_tint
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[7].default_value = mat.BDSF_1_roughness
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[8].default_value = mat.BDSF_1_anisotropic
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[9].default_value = mat.BDSF_1_anisotropic_rot
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[10].default_value = mat.BDSF_1_sheen
#   bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[11].default_value = mat.BDSF_1_sheen_tint
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[12].default_value = mat.BDSF_1_clearcoat
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[13].default_value = mat.BDSF_1_clear_rough
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[14].default_value = mat.BDSF_1_IOR
#   bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[15].default_value = mat.BDSF_1_transmission
#   bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[16].default_value = mat.BDSF_1_transmission_rough
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[17].default_value = mat.BDSF_1_emission_color
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[18].default_value = mat.BDSF_1_emission_strength
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF"].inputs[19].default_value = mat.BDSF_1_alpha
#   bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[0].default_value = mat.BDSF_2_hue
#    bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[1].default_value = mat.BDSF_2_saturation
#   bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[2].default_value = mat.BDSF_2_value
#    bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[3].default_value = mat.BDSF_2_fac
#    bpy.data.materials[mat.name].node_tree.nodes["Hue Saturation Value"].inputs[4].default_value = mat.BDSF_2_color
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[1].default_value = mat.BDSF_2_subsurface
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[2].default_value[0] = mat.BDSF_2_subsurf_radius[0]
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[2].default_value[1] = mat.BDSF_2_subsurf_radius[1]
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[2].default_value[2] = mat.BDSF_2_subsurf_radius[2]
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[3].default_value = mat.BDSF_2_subsurf_color
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[4].default_value = mat.BDSF_2_metallic
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[5].default_value = mat.BDSF_1_specular
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[6].default_value = mat.BDSF_2_specular_tint
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[7].default_value = mat.BDSF_2_roughness
 #   bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[8].default_value = mat.BDSF_2_anisotropic
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[9].default_value = mat.BDSF_2_anisotropic_rot
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[10].default_value = mat.BDSF_2_sheen
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[11].default_value = mat.BDSF_2_sheen_tint
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[12].default_value = mat.BDSF_2_clearcoat
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[13].default_value = mat.BDSF_2_clear_rough
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[14].default_value = mat.BDSF_2_IOR
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[15].default_value = mat.BDSF_2_transmission
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[16].default_value = mat.BDSF_2_transmission_rough
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[17].default_value = mat.BDSF_2_emission_color
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[18].default_value = mat.BDSF_2_emission_strength
#    bpy.data.materials[mat.name].node_tree.nodes["Principled BSDF.001"].inputs[19].default_value = mat.BDSF_2_alpha
#    bpy.context.object.active_material.blend_method = mat.blend_method
#    bpy.context.object.active_material.shadow_method = mat.shadow_method
#    bpy.context.object.active_material.refraction_depth = mat.refraction_depth
##    bpy.context.object.active_material.pass_index = mat.pass_index
#    bpy.context.object.active_material.diffuse_color = mat.diffuse_color
#    bpy.context.object.active_material.metallic = mat.metallic
#    bpy.context.object.active_material.roughness = mat.roughness
#    bpy.ops.node.add_node(type="ShaderNodeTexNoise", use_transform=True)
#    bpy.ops.node.select(wait_to_deselect_others=True, mouse_x=114, mouse_y=765, extend=False, deselect_all=True)
#    bpy.ops.node.link(detach=False)

def add_material_cell(mat_name, r, g, b):
    #adds a soap bubble like material for use in rendering cells
    print ("add cell material " + mat_name)
    
    # check whether the material already exists
    if bpy.data.materials.get(mat_name):
        mat = bpy.data.materials[mat_name]
    else:
        # create the material
        mat = bpy.data.materials.new(mat_name)
        
    mat.diffuse_color = (1.0, 1.0, 1.0, 1.0) # viewport color
    mat.use_nodes = True
    mat.blend_method = 'BLEND'

    # get the material nodes
    nodes = mat.node_tree.nodes

    # clear all nodes to start clean
    for node in nodes:
        nodes.remove(node)

    # create principled node for main color
    node_main = nodes.new(type='ShaderNodeBsdfPrincipled')
    node_main.location = -200,100
    node_main.inputs['Base Color'].default_value = (r,g,b,1)
    node_main.inputs['Metallic'].default_value = 0.136
    node_main.inputs['Specular'].default_value = 0.500
    node_main.inputs['Specular Tint'].default_value = 0.555
    node_main.inputs['Roughness'].default_value = 0.318
    node_main.inputs['Anisotropic'].default_value = 0.041
    node_main.inputs['Anisotropic Rotation'].default_value = 0.048
    node_main.inputs['Sheen'].default_value = 0.052
    node_main.inputs['Sheen Tint'].default_value = 0.030
    node_main.inputs['Clearcoat'].default_value = 0.114
    node_main.inputs['Clearcoat Roughness'].default_value = 0.123
    node_main.inputs['IOR'].default_value = 1.450
    node_main.inputs['Transmission'].default_value = 0.882
    node_main.inputs['Transmission Roughness'].default_value = 0.0
    node_main.inputs['Alpha'].default_value = 0.414
    
    # create noise texture source
    node_noise = nodes.new(type="ShaderNodeTexNoise")
    node_noise.inputs['Scale'].default_value = 0.600
    node_noise.inputs['Detail'].default_value = 15.0
    node_noise.inputs['Roughness'].default_value = 0.500
    node_noise.inputs['Distortion'].default_value = 3.0
    
    # create HSV 
    node_HSV = nodes.new(type="ShaderNodeHueSaturation")
    node_HSV.inputs['Hue'].default_value = 0.800
    node_HSV.inputs['Saturation'].default_value = 2.00
    node_HSV.inputs['Value'].default_value = 2.00
    node_HSV.inputs['Fac'].default_value = 1.00
   
    # create second principled node for random color variation
    node_random = nodes.new(type='ShaderNodeBsdfPrincipled')
    node_random.location = -200,-100
    node_random.inputs['Base Color'].default_value = (r,g,b,1) #Hue Saturation Value
    node_random.inputs['Metallic'].default_value = 0.0
    node_random.inputs['Specular'].default_value = 0.500
    node_random.inputs['Specular Tint'].default_value = 0.0
    node_random.inputs['Roughness'].default_value = 0.482
    node_random.inputs['Anisotropic'].default_value = 0.0
    node_random.inputs['Anisotropic Rotation'].default_value = 0.0
    node_random.inputs['Sheen'].default_value = 0.0
    node_random.inputs['Sheen Tint'].default_value = 0.0
    node_random.inputs['Clearcoat'].default_value = 0.0
    node_random.inputs['Clearcoat Roughness'].default_value = 0.0
    node_random.inputs['IOR'].default_value = 1.450
    node_random.inputs['Transmission'].default_value = 1.0
    node_random.inputs['Transmission Roughness'].default_value = 0.0
    node_random.inputs['Alpha'].default_value = 0.555

    # create mix shader node
    node_mix = nodes.new(type='ShaderNodeMixShader')
    node_mix.location = 0,0
    node_mix.inputs['Fac'].default_value = 0.079

    # create output node
    node_output = nodes.new(type='ShaderNodeOutputMaterial')   
    node_output.location = 200,0

    # link nodes
    links = mat.node_tree.links    
    link_noise_HSV = links.new(node_noise.outputs[1], node_HSV.inputs[4])
    link_HSV_random = links.new(node_HSV.outputs[0], node_random.inputs[0])
    link_main_mix = links.new(node_main.outputs[0], node_mix.inputs[1])
    link_random_mix = links.new(node_random.outputs[0], node_mix.inputs[2])
    link_mix_out = links.new(node_mix.outputs[0], node_output.inputs[0])
        

class Material():
    def __init__(self):
        self.name = "Cell Material"
        self.BDSF_1_fac = 0.079
        self.BDSF_1_color = (0.00749999, 0.020955, 0.3, 1)
        self.BDSF_1_subsurface = 0
        self.BDSF_1_subsurf_radius = (1.1, 0.2, 0.2)
        self.BDSF_1_subsurf_color = (0.8, 0.8, 0.8, 1)
        self.BDSF_1_metallic = 0.136
        self.BDSF_1_specular = 0.5
        self.BDSF_1_specular_tint = 0.555
        self.BDSF_1_roughness = 0.318
        self.BDSF_1_anisotropic = 0.041
        self.BDSF_1_anisotropic_rot = 0.048
        self.BDSF_1_sheen = 0.052
        self.BDSF_1_sheen_tint = 0.03
        self.BDSF_1_clearcoat = 0.114
        self.BDSF_1_clear_rough = 0.123
        self.BDSF_1_IOR = 1.45
        self.BDSF_1_transmission = 0.882
        self.BDSF_1_transmission_rough = 0
        self.BDSF_1_emission_color = (0, 0, 0, 1)
        self.BDSF_1_emission_strength = 1
        self.BDSF_1_alpha = 0.414

        self.BDSF_2_hue = 0.8
        self.BDSF_2_saturation = 2
        self.BDSF_2_value = 2
        self.BDSF_2_fac = 1
        self.BDSF_2_color = (0.8, 0.8, 0.8, 1)
        self.BDSF_2_subsurface = 0
        self.BDSF_2_subsurf_radius = (1.1, 0.2, 0.2)
        self.BDSF_2_subsurf_color = (0.8, 0.8, 0.8, 1)
        self.BDSF_2_metallic = 0
        self.BDSF_2_specular = 0.5
        self.BDSF_2_specular_tint = 0
        self.BDSF_2_roughness = 0.482
        self.BDSF_2_anisotropic = 0
        self.BDSF_2_anisotropic_rot = 0
        self.BDSF_2_sheen = 0
        self.BDSF_2_sheen_tint = 0
        self.BDSF_2_clearcoat = 0
        self.BDSF_2_clear_rough = 0
        self.BDSF_2_IOR = 1.45
        self.BDSF_2_transmission = 1
        self.BDSF_2_transmission_rough = 0
        self.BDSF_2_emission_color = (0, 0, 0, 1)
        self.BDSF_2_emission_strength = 1
        self.BDSF_2_alpha = 0.555

        self.blend_method = 'BLEND'
        self.shadow_method = 'NONE'
        self.refraction_depth = 0
        self.pass_index = 0
        self.diffuse_color = (0.8, 0.8, 0.8, 1)
        self.metallic = 0
        self.roughness = 0.4

def initialize_cells(num_cells, loc_array, material):
    if len(num_cells) != len(loc_array):
        print("Number of cells must match number of cell locations")
        return
    add_material(material)
    for i in range(num_cells):
        cell = Cell("cell_" + str(i), loc = loc_array[i])
        make_cell(cell)

def setup_world():
    print("setup worldDDDD")
    bpy.context.scene.use_gravity = False
    bpy.context.scene.unit_settings.system = 'METRIC'
    bpy.context.scene.unit_settings.scale_length = 1
    bpy.context.scene.unit_settings.system_rotation = 'DEGREES'
    bpy.context.scene.unit_settings.length_unit = 'METERS'
    bpy.context.scene.unit_settings.mass_unit = 'KILOGRAMS'
    bpy.context.scene.unit_settings.time_unit = 'SECONDS'
    bpy.context.scene.unit_settings.temperature_unit = 'KELVIN'

    #add materials
    
    add_world_HDRI()
  
def add_world_HDRI():
    print("add world HDRI")
    C = bpy.context
    scn = C.scene

    # Get the environment node tree of the current scene
    node_tree = scn.world.node_tree
    tree_nodes = node_tree.nodes

    # Clear all nodes
    tree_nodes.clear()

    # Add Background node
    node_background = tree_nodes.new(type='ShaderNodeBackground')

    # Add Environment Texture node
    node_environment = tree_nodes.new('ShaderNodeTexEnvironment')
    # Load and assign the image to the node property
    node_environment.image = bpy.data.images.load("//missile_launch_facility_01_4k.hdr") # Relative path- must be in same directory as blend file
    node_environment.location = -300,0

    # Add Output node
    node_output = tree_nodes.new(type='ShaderNodeOutputWorld')   
    node_output.location = 200,0

    # Link all nodes
    links = node_tree.links
    link = links.new(node_environment.outputs["Color"], node_background.inputs["Color"])
    link = links.new(node_background.outputs["Background"], node_output.inputs["Surface"])
    
    #set film to transparent to hide background
    bpy.context.scene.render.film_transparent = True
    
    #change render preview mode
    for area in bpy.context.screen.areas:   #only updates windows in current tab, e.g. Sxripting but not Layout
        if area.type == 'VIEW_3D':
            print ("update view 3d to rendered")
            space = area.spaces.active
            if space.type == 'VIEW_3D':
                space.shading.type = 'RENDERED'
        
def render(file_path,scene,start,end):
    scene.render.image_settings.file_format = 'PNG'
    old_fp = scene.render.filepath
    scene.render.filepath = file_path
    scene = bpy.context.scene
    scene.frame_start = start
    scene.frame_end = end
    handlers = bpy.app.handlers.frame_change_post.copy()
    bpy.app.handlers.frame_change_post.clear()
    for frame in range(scene.frame_start, scene.frame_end):
        bpy.context.scene.frame_set(frame)
        for func in handlers:
            func(scene)
        file_name = "frame" + str(scene.frame_current)
        scene.render.filepath += file_name
        bpy.ops.render.render(write_still=True) 
        scene.render.filepath = scene.render.filepath.removesuffix(file_name)
    scene.render.filepath = old_fp
    for func in handlers:
        bpy.app.handlers.frame_change_post.append(func)
    return

class handler_class:
    def __init__(self):
        self.cell_types = ["Sphere","type1","type2"]
        self.division_rates = {}
        self.growth_rates = {}
        self.adhesion_forces = {} # dictionary of dictionaries 
                                  #ex: {"sphere": {"sphere": 100, "type1": 100, "type2": 100},"type1": {"sphere": 200, "type1": 200, "type2": 200},"type2": {"sphere": 300, "type1": 300, "type2": 300}}
                                  # in the example above, spheres exert 100 force on other cells, type1 exerts 200, type2 exerts 300
                                  # could change values to have different amounts of force on different cell types
        for type in self.cell_types:
            self.division_rates[type] = 0
            self.growth_rates[type] = 0
            self.adhesion_forces[type] = {}
            for i in self.cell_types:
                self.adhesion_forces[type][i] = 0
        self.active_cell_types = [] # add active types to know what collections to divide
        return
    def set_division_rate(self,cell_type,rate):
        # assume 60 frames per second 
        # rate is in divisions per second
        self.division_rates[cell_type] = 60/rate
        return
    def set_growth_rate(self,cell_type,rate):
        self.growth_rates[cell_type] = rate
        return
    def set_adhesion(self,type1, force_type1_to_type2, type2, force_type2_to_type1):
        self.adhesion_forces[type1][type2] = force_type1_to_type2
        self.adhesion_forces[type2][type1] = force_type2_to_type1
        return
    def div_handler(self,scene,depsgraph):
        for cell_type in self.active_cell_types:
            num_cells = len(bpy.data.collections[cell_type].objects)
            current_frame = bpy.data.scenes[0].frame_current
            if current_frame % self.division_rates[cell_type] == 0:
                print("DIVIDING CELLS")
                for i in range(num_cells):
                    cell_name = bpy.data.collections[cell_type].objects[i].name
                    cell = bpy.data.objects[cell_name]
                    bpy.data.objects[cell_name].select_set(True)
                    bpy.context.view_layer.objects.active = bpy.data.objects[cell_name]
                    print(cell_name)
                    turn_off_physics()
                    d1, d2 = divide(cell)
                    bpy.data.objects[d1.data["name"]].select_set(True)
                    bpy.data.objects[d2.data["name"]].select_set(False)
                    bpy.context.view_layer.objects.active = bpy.data.objects[d1.data["name"]]
                    turn_on_physics()
                    bpy.data.objects[d1.data["name"]].select_set(False)
                    bpy.data.objects[d2.data["name"]].select_set(True)
                    bpy.context.view_layer.objects.active = bpy.data.objects[d2.data["name"]]
                    turn_on_physics()
                    bpy.data.objects[d2.data["name"]].select_set(False)
    def growth_handler(self,scene, depsgraph): # WIP
        print("Frame:",scene.frame_current)
        cell_objs = [obj for obj in scene.objects if obj.name.startswith("Cell_")] # change between scene and depsgragh here
        
        num_cells = len(bpy.data.collections["Cells"].objects)
        for cell_obj in cell_objs:
            #cell_name = bpy.data.collections["Cells"].objects[i].name
            cell_name = cell_obj.name
            #cell_obj = bpy.data.objects[cell_name]
            #bpy.data.objects[cell_name].select_set(True)
            #bpy.context.view_layer.objects.active = bpy.data.objects[cell_name]
            #bpy.ops.object.modifier_apply(modifier="CLOTH")
            
            print(cell_obj.name," Volume:",Goo.calculate_volume(cell_obj)," Shrinking Factor:",cell_obj.modifiers["Cloth"].settings.shrink_min)
            cell_obj.modifiers["Cloth"].settings.shrink_min -= 0.01 # constantly changing shrink_min
