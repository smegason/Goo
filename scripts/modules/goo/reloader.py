import sys


def reset():
    to_delete = []
    for modname, _ in sys.modules.items():
        if modname.startswith("goo"):
            to_delete.append(modname)
    for modname in to_delete:
        del sys.modules[modname]
