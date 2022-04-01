# add_cell_controller

# from numpy import size

class add_cell_controller:
    def __init__(self):
        self.cell_types = ["sphere", "type1", "type2"]
        self.active_types = []
        self.active_cells = []
        return
    def remove_cell(self,index):
        t = self.active_cells[index].type
        self.active_cells.remove(self.active_cells[index])
        found = False
        for i in self.active_cells:
            if i.type == t:
                found = True
                break
        if found:
            return False, None
        row = self.active_types.index(t)
        self.active_types.remove(t)
        return True, row

    def add_cell(self, size, location, type):
        c = cell(size, location, type)
        self.active_cells.append(c)
        return c

    def check_location():
        # return false if the location would collide with another cell
        #       true otherwise
        return


class cell:
    def __init__(self, size, location, type):
        self.size = size
        self.location = location
        self.type = type
        self.name = type + " " + str(location)
        return
    