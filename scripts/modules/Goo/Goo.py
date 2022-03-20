# Goo.py - This is the Goo library. It contains all helper functions for Goo
# Goo is licensed under BSDv2

import bpy, bmesh, mathutils, numpy as np

# This function calculates the volume of a cell 
def calculate_volume(obj): # obj = bpy.context.object
    # We need to get the cell as it is evaluated in the simulation.
    # To do this, we fetch its dependency graph and obtain the
    # evaluated cell (denoted as obj_eval here)
    dg = bpy.context.evaluated_depsgraph_get()
    obj_eval = obj.evaluated_get(dg)
    # We use this graph to obtain a new mesh from which
    # we can calculate the volume using Blender's bmesh library
    mesh_from_eval = obj_eval.to_mesh()
    bm = bmesh.new()
    bm.from_mesh(mesh_from_eval)
    volume = bm.calc_volume(signed = True)
    return volume

# In standard Goo simulations, cells divide along the major axis.
# This function finds the major axis of a cell using Principle
# Component Analysis. 
def get_major_axis(obj): # obj = bpy.context.object
    # We need to get the cell as it is evaluated in the simulation.
    # To do this, we fetch its dependency graph and obtain the
    # evaluated cell (denoted as obj_eval here)
    dg = bpy.context.evaluated_depsgraph_get()
    obj_eval = obj.evaluated_get(dg)
    # We obtain the (x,y,z) coordinates of the vertices in the
    # evaluated cell
    vertices = obj_eval.data.vertices
    vert_coords = [(obj_eval.matrix_world @ v.co) for v in vertices]
    vert_coords = np.asarray(vert_coords)

    # We separate the x, y, and z coordinates into their own arrays.
    # We also subtract the mean of each dimension from the corresponding
    # array values. This is part of the PCA algorithm.
    x = vert_coords[:, 0]
    x = x - np.mean(x)
    y = vert_coords[:, 1]
    y = y - np.mean(y)
    z = vert_coords[:, 2]
    z = z - np.mean(z)

    # We stack the three arrays together to make the "new" coordinates
    new_coords = np.vstack([x, y, z])
    # This is then used to find the covariance matrix of the coordinates
    cov_matrix = np.cov(new_coords)
    # Per the PCA algorithm, we find the eigenalues and eigenvectors 
    # of the covariance matrix
    eigenvals, eigenvecs = np.linalg.eig(cov_matrix)

    # The eigenvalues are sorted, and the primary eigenvector
    # is the major axis.
    sort_indices = np.argsort(eigenvals)
    major_x, major_y, major_z = eigenvecs[:, sort_indices[-1]]
    major_axis = (major_x, major_y, major_z)
    return major_axis

# This function returns 2 angles associated with the
# division axis of a cell: The angle angle between the the division axis
# and the z-axis (phi, in spherical coordinates) and the angle
# projected on the xy-plan (theta)
def get_division_angles(axis):
    # We define the unit vector of z-axis 
    z_axis = np.array((0, 0, 1))
    # We find the unit vector of the division axis
    division_axis = axis/np.linalg.norm(axis)
    # We calculate the dot product of the normalized
    # z and division axes
    dot_product = np.dot(z_axis, division_axis)
    # The inverse cosine of the dot product is phi
    phi = np.arccos(dot_product)

    # The first step to find theta is to find the projection
    # of the division axis on the z axis
    proj_axis_on_z = dot_product*z_axis
    # This projection is subtracted from the division axis
    # To find the projection on the xy-plane
    proj_xy_axis = division_axis - proj_axis_on_z
    # We normalize this projection
    proj_xy_axis = proj_xy_axis/np.linalg.norm(proj_xy_axis)
    # We take the dot product of the x-axis and the normalized projection
    dot_product = np.dot((1, 0, 0), proj_xy_axis)
    # The inverse cosin of this dot product is theta.
    theta = np.arccos(dot_product)

    #prev_axis_norm = prev_axis/np.linalg.norm(prev_axis)
    #dot_product = np.dot(division_axis, prev_axis_norm)
    #div_angle_change = np.arccos(dot_product)
    #return phi, theta, div_angle_change
    print("GOT DIVISION ANGLES")
    return phi, theta

# This function calculates the contact area of a cell
# using the angles between adjacent faces
def calculate_contact_area(obj): # obj = bpy.context.object
    # We need to get the cell as it is evaluated in the simulation.
    # To do this, we fetch its dependency graph and obtain the
    # evaluated cell (denoted as obj_eval here)
    dg = bpy.context.evaluated_depsgraph_get()
    obj_eval = obj.evaluated_get(dg)
    # We use this graph to obtain a new mesh from which
    # we can calculate the angles between adjacent faces of the cell
    mesh_from_eval = obj_eval.to_mesh()
    bm = bmesh.new()
    bm.from_mesh(mesh_from_eval)

    # The contact area variable is initialized to zero.
    contact_area = 0
    # All faces on the mesh are initially deselected
    for face in bm.faces:
        face.select = False

    # We loop over all the edges in the mesh
    for edge in bm.edges:
        # We find the angle between the two adjacent phases
        # of the edge. If the angle is less than 0.01, the
        # adjacent phases are either flat or concave, so this
        # face is presumed to be in contact with another object 
        # and is selected
        angle = edge.calc_face_angle_signed()
        if angle < 0.01:
            for face in edge.link_faces:
                face.select = True

    # We loop through the selected faces and add the area of
    # each face to the contact area variable            
    for face in bm.faces:
        if face.select == True:
            contact_area += face.calc_area()

    return contact_area

# The repair hole function adds a new face to the cell after
# division of the mesh and regularizes the mesh so that the
# mesh is evenly covered with faces
def repair_hole(obj): 
    #bpy.context.view_layer.objects.active = obj
    # Go into Blender Edit Mode
    bpy.ops.object.mode_set(mode = 'EDIT')
    # Select all the edges in the mesh
    bpy.ops.mesh.select_mode(type="EDGE") 
    bpy.ops.mesh.select_all(action = 'SELECT')
    # Add a new face (based on open eduges)
    bpy.ops.mesh.edge_face_add()
    # Deselct the edges
    bpy.ops.mesh.select_all(action = 'DESELECT')
    # Return to Blender Object Mode
    bpy.ops.object.mode_set(mode = 'OBJECT')
    # Remesh the cell with voxels of size 0.2
    bpy.ops.object.modifier_add(type='REMESH')
    bpy.context.object.modifiers["Remesh"].voxel_size = 0.2
    bpy.ops.object.modifier_apply(modifier="Remesh")

# The sep function splits the cell mesh in two for division
def sep(obj):
    # We need to get the cell as it is evaluated in the simulation.
    # To do this, we fetch its dependency graph and obtain the
    # evaluated cell (denoted as obj_eval here)
    bpy.ops.object.mode_set(mode='OBJECT')
    dg = bpy.context.evaluated_depsgraph_get()
    obj = obj.evaluated_get(dg)
    # The mother cell name is the name of the cell currently being divided
    mother_name = obj.name
    # By Blender convention, duplicates of existing objects are given
    # the same name as the original object, but with ".001" added to the end
    # We use the naming convention to tentatively set the name of one 
    # daughter cell
    daughter_name = mother_name + ".001"
    # Get the cell's center of mass
    bpy.ops.object.origin_set(type='ORIGIN_CENTER_OF_MASS')
    COM = obj.location
    # Go into edit mode
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='SELECT')
    # Get the cell's major axis
    major_axis = get_major_axis(obj)
    # Get the division angles, if desired
    #phi, theta = get_division_angles(major_axis)
    #print(mother_name + " PHI: " + str(phi) + ", THETA: " + str(theta))
    # Bisect the cell along the plane specified by the center of mass (point)
    # and major axis (normal vectorr)
    bpy.ops.mesh.bisect(plane_co = COM, plane_no = major_axis, use_fill=False, flip=False)
    # Go into object mode
    bpy.ops.object.mode_set(mode = 'OBJECT')
    # Now we separate the split mesh into two separate objects.
    # First, obtain the vertex coordinates of the bisected mother cell.
    obj = bpy.data.objects[mother_name]
    new_vertices = obj.data.vertices
    new_vert_coords = [(obj.matrix_world @ v.co) for v in new_vertices]
    new_vert_coords = np.asarray(new_vert_coords)

    # We will choose half of the vertices to be separated into a new object.
    # We create a list of the vertex indices we will use to create the
    # daughter cell object
    separation_indices = []
    # We loop through each vertex and calculate the signed distance between the vertex
    # and the division plane (which is specified by the center of mass and major axis).
    # If the distance is greater than -0.05 (the vertex is on a specific half of the cell),
    # we add the vertex's index i to our list of separation indices
    for i in range(len(new_vert_coords)):
        distance = mathutils.geometry.distance_point_to_plane(new_vert_coords[i], COM, major_axis)
        if distance > -0.05:
            separation_indices.append(i)
    # We go back into edit mode and deselect all vertices
    bpy.ops.object.mode_set(mode = 'EDIT') 
    bpy.ops.mesh.select_mode(type="VERT")
    bpy.ops.mesh.select_all(action = 'DESELECT')
    # We return to object mode and loop through all the vertices.
    # If the vertex's index ix contained in the list of separation indices,
    # we select it for separation
    bpy.ops.object.mode_set(mode = 'OBJECT')
    for index in separation_indices:
        obj.data.vertices[index].select = True
    # We go back into edit mode and separate the selected vertices
    # as a new daughter cell.
    bpy.ops.object.mode_set(mode = 'EDIT') 
    bpy.ops.mesh.separate(type='SELECTED')
    # We return to object mode and select only the "original" mother cell
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
    
    #smooth mesh
    bpy.ops.object.select = True
    bpy.ops.object.shade_smooth()
    
    #add subsurface modifier to make smoother
    bpy.ops.object.modifier_add(type = 'SUBSURF')
    bpy.context.object.modifiers["Subdivision"].levels = cell.data['subdiv']
    
    #add cloth settings for physics
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
    
    #add Collision modifier for physics
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
    
    #add material to cell based on name of material
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

class Force():
    def __init__(self, force_name, cell_name, strength):
        self.name = force_name
        self.strength = strength
        self.associated_cell = cell_name
        self.falloff_power = 0
    def get_blender_force(self):
        obj = bpy.data.objects[self.name]
        return obj

def make_force(force):
    bpy.ops.object.effector_add(type='FORCE', enter_editmode=False, align='WORLD', location=bpy.data.objects[force.associated_cell].location, scale=(1, 1, 1))
    bpy.context.object.field.strength = force.strength
    bpy.context.object.name = force.name
    bpy.context.object.field.falloff_power = force.falloff_power

def initialize_cells(num_cells, loc_array, material):
    if len(num_cells) != len(loc_array):
        print("Number of cells must match number of cell locations")
        return
    add_material(material)
    for i in range(num_cells):
        cell = Cell("cell_" + str(i), loc = loc_array[i])
        make_cell(cell)

def setup_world():
    """
    This function sets up the default values used for simulations in Goo including units and rendering background
    
    :return: None
    """    
    print("setup world")
    bpy.context.scene.use_gravity = False
    bpy.context.scene.unit_settings.system = 'METRIC'
    bpy.context.scene.unit_settings.scale_length = 1
    bpy.context.scene.unit_settings.system_rotation = 'DEGREES'
    bpy.context.scene.unit_settings.length_unit = 'MICROMETERS'
    bpy.context.scene.unit_settings.mass_unit = 'MILLIGRAMS'
    bpy.context.scene.unit_settings.time_unit = 'SECONDS'
    bpy.context.scene.unit_settings.temperature_unit = 'CELSIUS'

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
    scripts_paths = bpy.utils.script_paths()
    node_environment.image = bpy.data.images.load(scripts_paths[-1]+"/modules/Goo/missile_launch_facility_01_4k.hdr") # Relative path- this file must be in same directory as blend file
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
