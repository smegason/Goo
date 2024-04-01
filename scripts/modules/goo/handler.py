# TODO: this is probably not necessary, handlers can exist within motion, divider submodules
class Handler:
    def run(self, scene, depsgraph):
        raise NotImplementedError("Subclasses must implement run() method.")


# TODO: handler to update mesh origins to COM
class DivisionHandler(Handler):
    pass


class MotionHandler(Handler):
    pass


class ForceUpdateHandler(Handler):
    # TODO: change to forces rather than cells. Adhesion forces will move with the cells themselves.
    def setup(self, get_cells, dt):
        self.get_cells = get_cells
        self.dt = dt

    # TODO: change force size to match cell size
    def run(self, scene, depsgraph):
        for cell in self.get_cells():
            for force in cell.forces:
                force.loc = cell.COM()


class TimingHandler(Handler):
    def __init__(self, start_time):
        self.start_time = start_time

    def run(self):
        pass
