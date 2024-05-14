import numpy as np
import re # Regular Expressions
import os
import utilities as deps
from pathlib import Path


# Change the working directory
filepath = r"D:\Actual Documents\UC Davis\Programming\VSP2AVL\test_data\PelicanC6_DegenGeom.csv" # CSV file name
savepath = r"D:\Actual Documents\UC Davis\Programming\VSP2AVL\save_data"
Sref = 8500
Cref = 33
Bref = 265
Xref, Yref, Zref = [122.338097-0.048*Cref, 0, 2.78]
mach_number = 0.82 #input("Input default mach number: ")
tolerance = 0.05
write_bodies = True
vortices_per_unit_length = 0.5

# change directory and specify file name
path = Path(filepath).absolute()
os.chdir(path.parent)
filename = path.name
print('Directory changed to "{}"'.format(path.parent))
print('Now translating "{}"'.format(filename))

# airfoil_location = 'data/'
AVL_filename = re.split(r'\.',filename)[0] # sets AVL filename to same as DegenGeom file
component_delimiter = '# DegenGeom Type, Name, SurfNdx' # delimiter between component section of CSV


with open(filename) as f:
    DegenGeom = f.readlines()

components = []

# Look through CSV to find the starting line of each component. Save starting line index if it is a lifting surface, the name, and the ID of the component
for i, line in enumerate(DegenGeom):
    if component_delimiter in line:
        component = {}
        component['begin_index'] = i

        next_line = re.split(r',\s*',DegenGeom[i+1]) # thank you Weston
        component['is_lifting_surface'] = next_line[0] == 'LIFTING_SURFACE'
        component['is_body'] = next_line[0] == 'BODY'

        component['name'] = next_line[1]
        component['num'] = next_line[2]
        component['ID'] = next_line[3]

        components.append(component)

# For each component, if it is a lifting surface, find the STICK_NODE section and save the beginning and end of that section
for i, component in enumerate(components):
    print(component['name'] + ' ' + str(component['num']))
    # set ending index for component
    if i < len(components) - 1: #  if not last component in file
        components[i]['end_index'] = components[i+1]['begin_index'] # set ending index as beginning index of next component
    else: # if last component in file
        components[i]['end_index'] = len(DegenGeom)


    if not component['is_lifting_surface']:
        components[i] = deps.get_SURFACE_NODE_data(component, DegenGeom)
    
    if component['is_lifting_surface'] or not component['standard_body']:
        components[i] = deps.get_STICK_NODE_data(component, DegenGeom, tolerance)

    if component['is_lifting_surface']:
        components[i] = deps.get_control_surface_data(component, DegenGeom)
        components[i] = deps.interpret_control_surface(component, tolerance)

# data at top of AVL file
preamble = '''{}

#Mach
 {:.2f}

#IYsym   IZsym   Zsym
 0       0       0.0

#Sref    Cref    Bref
{:.2f}   {:.2f}    {:.2f}

#Xref    Yref    Zref
{:.2f}     {:.2f}     {:.2f}

'''.format(AVL_filename, mach_number, Sref, Cref, Bref, Xref, Yref, Zref)

AVL_file = [preamble]

# lifting component preamble
for component in components:
    if component['is_lifting_surface'] == True or (component['standard_body'] == False and write_bodies):
        AVL_file += deps.create_lifting_surface(component, vortices_per_unit_length)

    if component['is_body'] == True and component['standard_body'] == True and write_bodies:
        AVL_file += deps.create_body(component)

# write AVL file array to new file
with open(AVL_filename + '.avl', 'w+') as f:
    for string in AVL_file:
        f.write(string)
        

print("[VSP2AVL] AVL geometry file saved as \"{}.avl\"".format(AVL_filename))