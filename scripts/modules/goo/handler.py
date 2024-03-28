# TODO: this is probably not necessary, handlers can exist within motion, divider submodules
class Handler:
    def run(self):
        raise NotImplementedError("Subclasses must implement run() method.")


class DivisionHandler(Handler):
    pass


class MotionHandler(Handler):
    pass


class TimingHandler(Handler):
    def __init__(self, start_time):
        self.start_time = start_time

    def run(self):
        pass
