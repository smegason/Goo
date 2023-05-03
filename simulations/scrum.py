
    def set_random_motion_speed(self, motion_speed: float):
        self.random_motion_speed = motion_speed
    '''
    def get_bounding_box(obj):
        """returns the corners of the bounding box of an object in world coordinates
        """
        return [obj.matrix_world @ Vector(corner) for corner in obj.bound_box]
        
    def objectsOverlap(self, obj1, obj2):
        """returns True if the object's bounding boxes are overlapping
        """
        vert1 = self.worldBoundingBox(obj1)
        vert2 = self.worldBoundingBox(obj2)
        faces = [(0, 1, 2, 3), (4, 7, 6, 5), (0, 4, 5, 1), (1, 5, 6, 2), (2, 6, 7, 3), (4, 0, 3, 7)]

        bvh1 = BVHTree.FromPolygons(vert1, faces)
        bvh2 = BVHTree.FromPolygons(vert2, faces)

        return bool(bvh1.overlap(bvh2))
    '''

    def calculate_com(cell): 
        """Function to calculate the center of mass of a cell Blender object. 

        :param bpy.data.objects[cell.name] cell: cell
        :returns: The coordinate of the center of mass of the cell. 
        :rtype: Tuple(x,y,z)
        """
        bpy.context.view_layer.objects.active = cell
        dg = bpy.context.evaluated_depsgraph_get()
        cell_eval = cell.evaluated_get(dg)
        vertices = cell_eval.data.vertices
        vert_coords = np.asarray([(cell_eval.matrix_world @ v.co) for v in vertices])

        x = vert_coords[:, 0]
        y = vert_coords[:, 1]
        z = vert_coords[:, 2]
        COM = (np.mean(x), np.mean(y), np.mean(z))
        return COM



    def data_handler(self, scene, depsgraph): 

        # initialization
        location_list = []
        com_list = []
        distances = set()
        total_dist = 0
        current_frame = bpy.data.scenes[0].frame_current

        print(f'Frame number: {current_frame}')
        print(self.data_file_path)

        for collection in bpy.data.collections:
            if 'Cells' in collection.name_full:
                coll_name = collection.name_full

                for idx, cell in enumerate(collection.objects):
                    location_list.append((cell.name, cell.location))
                    COM = self.calculate_com(cell)
                    
                    # initialize the key for the cell name if the class dictionary 
                    com_list.append(COM)

                for i in range(len(com_list)):
                    for j in range(len(com_list)):
                        # Calculate the Euclidean distance between the two coordinates
                        distance = math.sqrt(sum([(a - b) ** 2 for a, b in zip(com_list[i], com_list[j])]))
                        # Add the distance to the list of distances
                        distances.add(distance)

                total_dist = sum(distances)
                self.distances_tot.append(total_dist)
                self.frames.append(current_frame)

        with open(self.data_file_path, "w") as file1:
            #file1.write(f'Frame number: {current_frame}\n')  
            #file1.write(f'Cell locations: {location_list}\n')       
            #file1.write(f'Cell COM: {com_list}\n')     
            #file1.write(f'Total distance: {total_dist}\n')
            file1.write(f'Export distance data: {self.distances_tot}\n')     
            file1.write(f'Export frame data: {self.frames}\n')        
            file1.write(f'Export time data: {self.times}\n') 
    
        print(f'Location for the cells: {location_list}')
        print(f'COM for the cells: {com_list}')
        print(f'Total distance between cells: {total_dist}')


    def random_motion(self): 

        for collection in bpy.data.collections:
            if 'Cells' in collection.name_full:
                for cell in collection.objects:
                    cell.select_set(True)
                    x = random.uniform(-self.random_motion_speed, self.random_motion_speed)
                    y = random.uniform(-self.random_motion_speed, self.random_motion_speed)
                    z = random.uniform(-self.random_motion_speed, self.random_motion_speed)
                    translation_coord = (x,y,z)
                    # translate the object the given value towards its corresponding axis
                    bpy.ops.transform.translate(value=translation_coord)
                    cell.select_set(False)
    '''
    # handlers should have a least one argument even if they don't use it
    def random_motion_handler(self, scene, depsgraph):

        for collection in bpy.data.collections:
            if 'Cells' in collection.name_full:
                coll_name = collection.name_full
                for cell in collection.objects:
                    cell.select_set(True)
                    x = random.uniform(-self.random_motion_speed, self.random_motion_speed)
                    y = random.uniform(-self.random_motion_speed, self.random_motion_speed)
                    z = random.uniform(-self.random_motion_speed, self.random_motion_speed)
                    translation_coord = (x,y,z)
                    # translate the object the given value towards its corresponding axis
                    bpy.ops.transform.translate(value=translation_coord)
                    cell.select_set(False)
    '''

    def motion_handler(self, scene, depsgraph):

        for collection in bpy.data.collections:
            if 'Cells' in collection.name_full:
                coll_name = collection.name_full
                for cell1 in collection.objects:
                    for cell2 in collection.objects: 
                        # Calculate the distance between the two objects
                        distance = (cell1.location - cell2.location).length
            
                        # If the distance is smaller than 0.1, stop moving the objects
                        if distance < 0.1:
                            return 
                        self.random_motion()
                        