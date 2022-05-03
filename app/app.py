# makes a GUI for defining cells

import sys

# from cgitb import handler
import tkinter as tk
from tkinter import ttk


# from cv2 import add
# from sqlalchemy import ForeignKeyConstraint
import add_cell_controller

root = tk.Tk()
# root.withdraw()

# have all the cells call divide and check their size

addCellController = add_cell_controller.add_cell_controller()


class Lotfi(tk.Entry):
    def __init__(self, master=None, **kwargs):
        self.var = tk.StringVar()
        tk.Entry.__init__(self, master, textvariable=self.var, **kwargs)
        self.old_value = ''
        self.var.trace('w', self.check)
        self.get, self.set = self.var.get, self.var.set

    def check(self, *args):
        # if self.get().isdecimal():
        try:
            float(self.get())
            # the current value is only digits; allow this
            self.old_value = self.get()
        except (Exception):
            # there's non-digit characters in the input; reject this
            self.set(self.old_value)
            print(Exception)


class Intry(tk.Entry):
    def __init__(self, master=None, **kwargs):
        self.var = tk.StringVar()
        tk.Entry.__init__(self, master, textvariable=self.var, **kwargs)
        self.old_value = ''
        self.var.trace('w', self.check)
        self.get, self.set = self.var.get, self.var.set

    def check(self, *args):
        # if self.get().isdecimal():
        try:
            int(self.get())
            # the current value is only digits; allow this
            self.old_value = self.get()
        except Exception:
            # there's non-digit characters in the input; reject this
            print("Exception=", sys.exc_info())
            self.set(self.old_value)


root.title("goo")
# root.geometry("500x600+10+20")
root.resizable(True, True)

main_width = 700
section_height = 100
borderColor = "blue"
borderThickness = 3
sectionBackground = "gray"

# Main sections of the app
add_cell_section = tk.Frame(root,
                            width=main_width,
                            height=section_height,
                            bg=sectionBackground,
                            pady=3,
                            highlightbackground=borderColor,
                            highlightthickness=borderThickness)
cell_division_section = tk.Frame(root,
                                 width=main_width,
                                 height=section_height,
                                 bg=sectionBackground,
                                 pady=3,
                                 highlightbackground=borderColor,
                                 highlightthickness=borderThickness)
growth_rate_section = tk.Frame(root,
                               width=main_width,
                               height=section_height,
                               bg=sectionBackground,
                               pady=3,
                               highlightbackground=borderColor,
                               highlightthickness=borderThickness)
cell_adhesion_section = tk.Frame(root,
                                 width=main_width,
                                 height=section_height,
                                 bg=sectionBackground,
                                 pady=3,
                                 highlightbackground=borderColor,
                                 highlightthickness=borderThickness)
render_section = tk.Frame(root,
                          width=main_width,
                          height=section_height,
                          bg=sectionBackground,
                          pady=3,
                          highlightbackground=borderColor,
                          highlightthickness=borderThickness)
generate_script_section = tk.Frame(root,
                                   width=main_width,
                                   height=section_height,
                                   bg=sectionBackground,
                                   pady=3,
                                   highlightbackground=borderColor,
                                   highlightthickness=borderThickness)

# Place main sections of the app
root.grid_rowconfigure(1, weight=1)
root.grid_columnconfigure(0, weight=1)

add_cell_section.pack(side=tk.TOP, pady=3)
cell_division_section.pack(side=tk.TOP, pady=3)
growth_rate_section.pack(side=tk.TOP, pady=3)
cell_adhesion_section.pack(side=tk.TOP, pady=3)
render_section.pack(side=tk.TOP, pady=3)
generate_script_section.pack(side=tk.TOP, pady=3)

# Create widgets for "Initial Cells" Section
padx = 3

tk.Label(add_cell_section,
         text="Initial Cells").grid(column=0, row=0, padx=(5, 632), columnspan=9)
cell_type_dropdown = ttk.Combobox(add_cell_section, text="")

x_lable = tk.Label(add_cell_section, text="x:")
y_lable = tk.Label(add_cell_section, text="y:")
z_lable = tk.Label(add_cell_section, text="z:")
x_input = Lotfi(add_cell_section, width=5)
y_input = Lotfi(add_cell_section, width=5)
z_input = Lotfi(add_cell_section, width=5)
size_lable = tk.Label(add_cell_section, text="Size:")
cell_size_slider = tk.Scale(add_cell_section, orient=tk.HORIZONTAL, from_=0, to=100)
cell_size_unit = tk.Label(add_cell_section, text="Unit")
add_cell_button = tk.Button(add_cell_section, text="Add Cell")
cells_list = tk.Listbox(add_cell_section, height=4)
scrollbar = tk.ttk.Scrollbar(
    add_cell_section,
    orient='vertical',
    command=cells_list.yview,
)
cells_list['yscrollcommand'] = scrollbar.set
remove_cell_button = tk.Button(add_cell_section, text="Remove Cell")

# PLace widgets for "Initial Cells" Section
cell_type_dropdown.grid(column=1, row=1, padx=padx)
x_lable.grid(row=1, column=2, padx=(padx, 0))
x_input.grid(row=1, column=3, padx=(0, padx))
y_lable.grid(row=1, column=4, padx=(padx, 0))
y_input.grid(row=1, column=5, padx=(0, padx))
z_lable.grid(row=1, column=6, padx=(padx, 0))
z_input.grid(row=1, column=7, padx=(0, padx))
size_lable.grid(row=1, column=8, padx=(padx, 0))
cell_size_slider.grid(row=1, column=9)
cell_size_unit.grid(row=1, column=10, padx=(0, padx))
add_cell_button.grid(row=2, column=1)
cells_list.grid(row=2, column=2, columnspan=6)
scrollbar.grid(row=2, column=6)
remove_cell_button.grid(row=2, column=8, columnspan=3)

# Create widgets for "Cell Division" Section
tk.Label(cell_division_section,
         text="Cell Division").grid(column=0, row=0, padx=(5, 762), columnspan=4)
division_labels = []
division_sliders = []

# Place widgets for "Cell Division" Section

# Create widgets for "Cell Growth" Section
tk.Label(growth_rate_section,
         text="Cell Growth").grid(column=0, row=0, padx=(5, 766), columnspan=4)
growth_labels = []
growth_sliders = []

# PLace widgets for "Cell Growth" Section

# Create widgets for "Cell Adhesion" Section
tk.Label(cell_adhesion_section,
         text="Cell Adhesion").grid(column=0, row=0, padx=(5, 754), columnspan=4)
adhesion_labels = []
adhesion_sliders = []
adhesion_dropdowns = []
adhesion_slider_values = []

# Place widgets for "Cell Adhesion" Section

# Create widgets for "Render" Section
checkbox_var = tk.IntVar()
checkbox = tk.Checkbutton(render_section, text="Rendering", variable=checkbox_var)
# Label(render_section,text="Rendering").grid(column=1,row=0,padx=(5,776),columnspan=5)
FilePathLabel = tk.Label(render_section, text="FilePath:")
FilePathInput = tk.Entry(render_section, width=30)
NumFramesLabel = tk.Label(render_section, text="Number of Frames:")
NumFrames = Intry(render_section, width=5)

# Place widgets for "Render" Section
checkbox.grid(column=0, row=0, padx=(5, 758), columnspan=5)
FilePathLabel.grid(row=1, column=1)
FilePathInput.grid(row=1, column=2)
NumFramesLabel.grid(row=1, column=3)
NumFrames.grid(row=1, column=4)

# Create widgets for "Generate Script" Section
tk.Label(generate_script_section,
         text="Generation").grid(column=0, row=0, padx=(5, 770), columnspan=4)
gen_script_button = tk.Button(generate_script_section, text="Generate Script")
save_script_button = tk.Button(generate_script_section, text="Copy Script")
run_script_button = tk.Button(generate_script_section, text="Run Script")

# Place widgets for "Generate Script" Section
gen_script_button.grid(row=1, column=1, padx=65)
save_script_button.grid(row=1, column=2, padx=65)
run_script_button.grid(row=1, column=3, padx=65)


class Script:
    def __init__(self):
        self.script = ""
        self.cells = []
        self.imports = ["from goo import goo",
                        "import importlib",
                        "import bpy",
                        "",
                        "importlib.reload(goo)"]
        return

    def generateScript(self):
        self.script = ""
        self.script += "# Import Libraries and setup world\n"
        for line in self.imports:
            self.script += line + "\n"
        self.script += "goo.setup_world()\n"  # setup

        self.script += "\n# Add cells and cell collections\n"

        # clear handlers before adding new ones
        self.script += "bpy.app.handlers.frame_change_post.clear()\n"

        self.script += "master_coll = bpy.context.view_layer.layer_collection\n"
        for type in addCellController.active_types:
            self.script += "collection = "
            self.script += "bpy.context.blend_data.collections.new(name='"+type+"')\n"
            self.script += "bpy.context.collection.children.link(collection)\n"
        for cell in addCellController.active_cells:
            self.script += "bpy.context.view_layer.active_layer_collection = "
            self.script += "bpy.context.view_layer.layer_collection.children['"
            self.script += cell.type
            self.script += "']\n"
            location = ("(" +
                        str(cell.location[0]) +
                        "," +
                        str(cell.location[1]) +
                        "," +
                        str(cell.location[2]) +
                        ")")
            self.script += 'goo.make_cell(goo.Cell(name_string="'
            self.script += cell.type + location
            self.script += '",'
            self.script += '                       loc=' + location + '))\n'

        self.script += "\n# Add handlers for division, growth and adhesion\n"
        self.script += "handlers = goo.handler_class()\n"
        divide = False
        grow = False
        adhesion = False
        for t in addCellController.active_types:
            self.script += "handlers.active_cell_types += ['"+t+"']\n"
            div_rate = division_sliders[addCellController.active_types.index(t)].get()
            if div_rate:
                self.script += "handlers.set_division_rate('"+t+"',"+str(div_rate)+")\n"
                divide = True
            growth_rate = growth_sliders[addCellController.active_types.index(t)].get()
            if growth_rate:
                grow = True
                self.script += "handlers.set_growth_rate('"
                self.script += t
                self.script += "', " + str(growth_rate) + ")\n"
            for type in addCellController.active_types:
                row = addCellController.active_types.index(t)
                force = int(adhesion_slider_values[row][type])
                if force != 0:
                    self.script += "handlers.set_adhesion('"
                    self.script += t+"', '" + type + "', " + str(force)+")\n"
                    adhesion = True
        if divide:
            self.script += "bpy.app.handlers.frame_change_post."
            self.script += "append(handlers.div_handler)\n"
        if grow:
            self.script += "bpy.app.handlers.frame_change_post."
            self.script += "append(handlers.growth_handler)\n"
        if adhesion:
            # make collections for forcefields
            self.script += "goo.make_force_collections"
            self.script += "(master_coll,handlers.active_cell_types)\n"
            self.script += "handlers.apply_forces()\n"
            self.script += "bpy.app.handlers.frame_change_post."
            self.script += "append(handlers.adhesion_handler)\n"

        if checkbox_var.get() == 1:
            self.script += "\n# Render animation\n"
            self.script += "scene = bpy.context.scene\n"
            # self.script += "fp = '/Users/michaelmitschjr/Desktop/Python/data/'\n"
            self.script += "fp = " + FilePathInput.get() + "\n"
            self.script += "scene.render.image_settings.file_format = 'PNG'\n"
            self.script += "goo.render(fp,scene,1,60)\n"

        print(self.script)
        return

    def saveScript(self):
        root.clipboard_clear()
        root.clipboard_append(self.script)
        root.update()  # now it stays on the clipboard after the window is closed
        return

    def runScript(self):
        return


script = Script()


def adhesion_dropdown_changed(event):
    dropdown = event.widget
    row = dropdown.grid_info()['row']-1
    key = dropdown.get()
    adhesion_sliders[row].set(adhesion_slider_values[row][key])


def adhesion_slider_changed(event):
    slider = event.widget
    row = slider.grid_info()['row']-1
    value = slider.get()
    key = adhesion_dropdowns[row].get()
    adhesion_slider_values[row][key] = int(value)
    print(value)
    print(row)


def addCell(event):
    location = (x_input.get(), y_input.get(), z_input.get())
    size = cell_size_slider.get()
    cellType = addCellController.cell_types[cell_type_dropdown.current()]
    cell = addCellController.add_cell(size, location, cellType)
    cells_list.insert(len(addCellController.active_cells), cell.name)
    script.cells += [cell]
    if cell.type in addCellController.active_types:
        return
    addCellController.active_types.append(cell.type)
    num_cells = len(addCellController.active_cells)
    division_labels.append(tk.Label(cell_division_section, text=cell.type))
    division_sliders.append(tk.Scale(cell_division_section,
                                     orient=tk.HORIZONTAL,
                                     length=500,
                                     from_=0,
                                     to=100))
    division_labels[-1].grid(column=1, row=num_cells, padx=20)
    division_sliders[-1].grid(row=num_cells, column=2, padx=10)

    growth_labels.append(tk.Label(growth_rate_section, text=cell.type))
    growth_sliders.append(tk.Scale(growth_rate_section,
                                   orient=tk.HORIZONTAL,
                                   length=500,
                                   from_=0,
                                   to=100))
    growth_labels[-1].grid(column=1, row=num_cells, padx=20)
    growth_sliders[-1].grid(row=num_cells, column=2, padx=10)

    adhesion_labels.append(tk.Label(cell_adhesion_section, text=cell.type))
    adhesion_dropdowns.append(tk.ttk.Combobox(cell_adhesion_section,
                                              text="",
                                              name=cell.type))
    adhesion_sliders.append(tk.Scale(cell_adhesion_section,
                                     orient=tk.HORIZONTAL,
                                     length=300,
                                     from_=-800,
                                     to=800))
    adhesion_labels[-1].grid(column=1, row=num_cells, padx=20)
    adhesion_dropdowns[-1].grid(column=2, row=num_cells, padx=20)
    adhesion_sliders[-1].grid(row=num_cells, column=3, padx=10)
    adhesion_dropdowns[-1].bind("<<ComboboxSelected>>", adhesion_dropdown_changed)
    adhesion_sliders[-1].bind("<ButtonRelease>", adhesion_slider_changed)
    forces = {}
    for t in addCellController.active_types:
        forces[t] = 0
    adhesion_slider_values.append(forces)
    update_adhesion_dropdown()
    return


def removeCell(event):
    index = cells_list.curselection()[0]
    addCellController.active_cells[index].type
    cells_list.delete(index)
    script.cells.pop(index)
    isTypeGone, row = addCellController.remove_cell(index)
    if isTypeGone:
        division_labels[row].destroy()
        division_labels.remove(division_labels[row])
        division_sliders[row].destroy()
        division_sliders.remove(division_sliders[row])

        growth_labels[row].destroy()
        growth_labels.remove(growth_labels[row])
        growth_sliders[row].destroy()
        growth_sliders.remove(growth_sliders[row])

        adhesion_labels[row].destroy()
        adhesion_labels.remove(adhesion_labels[row])
        adhesion_sliders[row].destroy()
        adhesion_sliders.remove(adhesion_sliders[row])
        adhesion_dropdowns[row].destroy()
        adhesion_dropdowns.remove(adhesion_dropdowns[row])
        # adhesion_slider_values.remove(adhesion_slider_values[row])
        update_adhesion_dropdown()
    return


def update_adhesion_dropdown():
    for dropdown in adhesion_dropdowns:
        dropdown['values'] = addCellController.active_types
    for t in addCellController.active_types:
        for dict in adhesion_slider_values:
            if not (t in dict):
                dict[t] = 0


def generateScript(event):
    script.generateScript()
    return


def saveScript(event):
    script.saveScript()
    return


def runScript(event):
    script.runScript()
    return


# add bindings to interface
cell_type_dropdown['values'] = addCellController.cell_types
remove_cell_button.bind("<Button-1>", removeCell)
add_cell_button.bind("<Button-1>", addCell)

gen_script_button.bind("<Button-1>", generateScript)
save_script_button.bind("<Button-1>", saveScript)

root.mainloop()

# forces that affect one type of cell should all be in the same collection
