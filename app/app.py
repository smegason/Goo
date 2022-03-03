from tkinter import *
from tkinter import ttk

from cv2 import add
from sqlalchemy import ForeignKeyConstraint
from add_cell_controller import add_cell_controller, cell

root = Tk()
#root.withdraw()

# have all the cells call divide and check their size 

addCellController = add_cell_controller()

class Lotfi(Entry):
    def __init__(self, master=None, **kwargs):
        self.var = StringVar()
        Entry.__init__(self, master, textvariable=self.var, **kwargs)
        self.old_value = ''
        self.var.trace('w', self.check)
        self.get, self.set = self.var.get, self.var.set

    def check(self, *args):
        #if self.get().isdecimal(): 
        try:
            float(self.get())
            # the current value is only digits; allow this
            self.old_value = self.get()
        except:
            # there's non-digit characters in the input; reject this 
            self.set(self.old_value)
class Intry(Entry):
    def __init__(self, master=None, **kwargs):
        self.var = StringVar()
        Entry.__init__(self, master, textvariable=self.var, **kwargs)
        self.old_value = ''
        self.var.trace('w', self.check)
        self.get, self.set = self.var.get, self.var.set

    def check(self, *args):
        #if self.get().isdecimal(): 
        try:
            int(self.get())
            # the current value is only digits; allow this
            self.old_value = self.get()
        except:
            # there's non-digit characters in the input; reject this 
            self.set(self.old_value)

root.title("Goo")
#root.geometry("500x600+10+20")
root.resizable(True,True)

main_width = 700
section_height = 100
borderColor = "blue"
borderThickness = 3
sectionBackground = "gray"

# Main sections of the app
add_cell_section = Frame(root,width=main_width,height=section_height,bg=sectionBackground,pady=3,highlightbackground=borderColor,highlightthickness=borderThickness)
cell_division_section = Frame(root, width=main_width, height=section_height,bg=sectionBackground,pady=3,highlightbackground=borderColor,highlightthickness=borderThickness)
growth_rate_section = Frame(root, width=main_width, height=section_height,bg=sectionBackground,pady=3,highlightbackground=borderColor,highlightthickness=borderThickness)
cell_adhesion_section = Frame(root, width=main_width, height=section_height, bg=sectionBackground,pady=3,highlightbackground=borderColor,highlightthickness=borderThickness)
render_section = Frame(root, width=main_width, height=section_height, bg=sectionBackground,pady=3,highlightbackground=borderColor,highlightthickness=borderThickness)
generate_script_section = Frame(root, width=main_width, height=section_height,bg=sectionBackground,pady=3,highlightbackground=borderColor,highlightthickness=borderThickness)

#Place main sections of the app
root.grid_rowconfigure(1, weight=1)
root.grid_columnconfigure(0, weight=1)

add_cell_section.pack(side=TOP,pady=3)
cell_division_section.pack(side=TOP,pady=3)
growth_rate_section.pack(side=TOP,pady=3)
cell_adhesion_section.pack(side=TOP,pady=3)
render_section.pack(side=TOP,pady=3)
generate_script_section.pack(side=TOP,pady=3)

#Create widgets for "Initial Cells" Section
padx = 3

Label(add_cell_section,text="Initial Cells").grid(column=0,row=0,padx=(5,632),columnspan=9)
cell_type_dropdown = ttk.Combobox(add_cell_section,text="")
x_lable = Label(add_cell_section,text="x:")
y_lable = Label(add_cell_section,text="y:")
z_lable = Label(add_cell_section,text="z:")
x_input = Lotfi(add_cell_section,width=5)
y_input = Lotfi(add_cell_section,width=5)
z_input = Lotfi(add_cell_section,width=5)
size_lable = Label(add_cell_section,text="Size:")
cell_size_slider = Scale(add_cell_section,orient=HORIZONTAL,from_=0,to=100)
cell_size_unit = Label(add_cell_section,text="Unit")
add_cell_button = Button(add_cell_section,text="Add Cell")
cells_list = Listbox(add_cell_section,height=4)
scrollbar = ttk.Scrollbar(
    add_cell_section,
    orient='vertical',
    command=cells_list.yview,
)
cells_list['yscrollcommand'] = scrollbar.set
remove_cell_button = Button(add_cell_section,text="Remove Cell")

#PLace widgets for "Initial Cells" Section
cell_type_dropdown.grid(column=1,row=1,padx=padx)
x_lable.grid(row=1,column=2,padx=(padx,0))
x_input.grid(row=1,column=3,padx=(0,padx))
y_lable.grid(row=1,column=4,padx=(padx,0))
y_input.grid(row=1,column=5,padx=(0,padx))
z_lable.grid(row=1,column=6,padx=(padx,0))
z_input.grid(row=1,column=7,padx=(0,padx))
size_lable.grid(row=1,column=8,padx=(padx,0))
cell_size_slider.grid(row=1,column=9)
cell_size_unit.grid(row=1,column=10,padx=(0,padx))
add_cell_button.grid(row=2,column=1)
cells_list.grid(row=2,column=2,columnspan=6)
scrollbar.grid(row=2,column=6)
remove_cell_button.grid(row=2,column=8,columnspan=3)

#Create widgets for "Cell Division" Section
Label(cell_division_section,text="Cell Division").grid(column=0,row=0,padx=(5,762),columnspan=4)
division_labels = []
division_sliders = []

#Place widgets for "Cell Division" Section

#Create widgets for "Cell Growth" Section
Label(growth_rate_section,text="Cell Growth").grid(column=0,row=0,padx=(5,766),columnspan=4)
growth_labels = []
growth_sliders = []

#PLace widgets for "Cell Growth" Section


#Create widgets for "Cell Adhesion" Section
Label(cell_adhesion_section,text="Cell Adhesion").grid(column=0,row=0,padx=(5,754),columnspan=4)
adhesion_labels = []
adhesion_sliders = []

#Place widgets for "Cell Adhesion" Section

#Create widgets for "Render" Section
checkbox_var = IntVar()
checkbox = Checkbutton(render_section,text="Rendering",variable=checkbox_var)
#Label(render_section,text="Rendering").grid(column=1,row=0,padx=(5,776),columnspan=5)
FilePathLabel = Label(render_section,text="FilePath:")
FilePathInput = Entry(render_section,width=30)
NumFramesLabel = Label(render_section,text="Number of Frames:")
NumFrames = Intry(render_section,width=5)

#Place widgets for "Render" Section
checkbox.grid(column=0,row=0,padx=(5,758),columnspan=5)
FilePathLabel.grid(row=1,column=1)
FilePathInput.grid(row=1,column=2)
NumFramesLabel.grid(row=1,column=3)
NumFrames.grid(row=1,column=4)

#Create widgets for "Generate Script" Section
Label(generate_script_section,text="Generation").grid(column=0,row=0,padx=(5,770),columnspan=4)
gen_script_button = Button(generate_script_section,text="Generate Script")
save_script_button = Button(generate_script_section,text="Copy Script")
run_script_button = Button(generate_script_section,text="Run Script")

#Place widgets for "Generate Script" Section
gen_script_button.grid(row=1,column=1, padx=65)
save_script_button.grid(row=1,column=2,padx=65)
run_script_button.grid(row=1,column=3,padx=65)

class Script:
    def __init__(self):
        self.script = ""
        self.cells = []
        self.imports = ["from Goo import Goo","import importlib","importlib.reload(Goo)","import bpy"]
        return
    def generateScript(self):
        self.script = ""
        for line in self.imports:
            self.script += line + "\n"
        self.script += "Goo.setup_world()\n" #setup
        self.script += "collection = bpy.context.blend_data.collections.new(name='Cells')\n" #add collection to store cells
        self.script += "bpy.context.collection.children.link(collection)\n"
        self.script += "layer_collection = bpy.context.view_layer.layer_collection.children[collection.name]\n" 
        self.script += "bpy.context.view_layer.active_layer_collection = layer_collection\n" # make the cells collection active
        self.script += "bpy.app.handlers.frame_change_post.clear()\n" #clear handlers before adding new ones
        for cell in self.cells:
            location = "("+str(cell.location[0])+","+str(cell.location[1])+","+str(cell.location[2])+")"
            self.script += 'cell = Goo.Cell(name_string = "cell_'+location+'", loc = '+location+')\n'
            self.script += "Goo.make_cell(cell)\n"
            #self.script += 'obj = bpy.data.objects["cell_"'+location+']]\n'
        if division_sliders[0].get():
            self.script += "bpy.app.handlers.frame_change_post.append(Goo.div_handler)\n"
        if checkbox_var.get() == 1:
            self.script += "scene = bpy.context.scene\n"
            #self.script += "fp = '/Users/michaelmitschjr/Desktop/Python/data/'\n"
            self.script += "fp = " + FilePathInput.get() + "\n"
            self.script += "scene.render.image_settings.file_format = 'PNG'\n"
            self.script += "Goo.render(fp,scene,1,60)\n"
        print(self.script)
        return

    def saveScript(self):
        root.clipboard_clear()
        root.clipboard_append(self.script)
        root.update() # now it stays on the clipboard after the window is closed
        return
    def runScript(self):

        return

script = Script()

def addCell(event):
    location = (x_input.get(),y_input.get(),z_input.get())
    size = cell_size_slider.get()
    cellType = addCellController.cell_types[cell_type_dropdown.current()]
    cell = addCellController.add_cell(size,location,cellType)
    cells_list.insert(len(addCellController.active_cells),cell.name)
    script.cells += [cell]
    if cell.type in addCellController.active_types:
        return
    addCellController.active_types.append(cell.type)
    num_cells = len(addCellController.active_cells)
    division_labels.append(Label(cell_division_section,text=cell.type))
    division_sliders.append(Scale(cell_division_section,orient=HORIZONTAL,length=500,from_=0,to=100))
    division_labels[-1].grid(column=1,row=num_cells,padx=20)
    division_sliders[-1].grid(row=num_cells,column=2,padx=10)

    growth_labels.append(Label(growth_rate_section,text=cell.type))
    growth_sliders.append(Scale(growth_rate_section,orient=HORIZONTAL,length=500,from_=0,to=100))
    growth_labels[-1].grid(column=1,row=num_cells,padx=20)
    growth_sliders[-1].grid(row=num_cells,column=2,padx=10)

    adhesion_labels.append(Label(cell_adhesion_section,text=cell.type))
    adhesion_sliders.append(Scale(cell_adhesion_section,orient=HORIZONTAL,length=500,from_=0,to=100))
    adhesion_labels[-1].grid(column=1,row=num_cells,padx=20)
    adhesion_sliders[-1].grid(row=num_cells,column=2,padx=10)
    return

def removeCell(event):
    index = cells_list.curselection()[0]
    t = addCellController.active_cells[index].type
    cells_list.delete(index)
    script.cells.pop(index)
    isTypeGone,row = addCellController.remove_cell(index)
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
    return

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
remove_cell_button.bind("<Button-1>",removeCell)
add_cell_button.bind("<Button-1>",addCell)

gen_script_button.bind("<Button-1>",generateScript)
save_script_button.bind("<Button-1>",saveScript)

root.mainloop()