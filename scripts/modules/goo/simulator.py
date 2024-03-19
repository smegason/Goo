import bpy
from handler import TimingHandler
import datetime


class Simulator:
    def __init__(self):
        self.handlers = []

    def add_handlers(self, handlers):
        self.handlers.extend(handlers)

    def add_handler(self, handler):
        self.handlers.append(handler)

    def run_simulation(self):
        start_time = datetime.now()
        bpy.app.handlers.frame_change_post.clear()
        bpy.app.handlers.frame_change_post.extend(self.handlers)
        bpy.app.handlers.frame_change_post.append(TimingHandler(start_time))
