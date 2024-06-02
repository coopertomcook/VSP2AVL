import re
import os
import numpy as np
import geom_data
from pathlib import Path
import matplotlib.pyplot as plt

#### USER SETTINGS #####################################
# Change the working directory
filepath = r"D:\Actual Documents\UC Davis\3rd Year\9 Spring 2024\EAE 130B\VSP2AVL\data\PelicanC6_DegenGeom.csv"                  # CSV file name
refNum = False                  # reference surface number (set to false if unknown)
Sref = 0                        # reference wing area
Cref = 0                        # reference chord length
Bref = 0                        # reference span
Xref, Yref, Zref = [0, 0, 0]    # center of gravity location
mach_number = 0.80              # default mach number
tolerance = 0.05                # minimum allowed geometric distance between sections for hingeline section creation (make small)
control_surfaces = True         # check for control surfaces in DegenGeom and add them to AVL file (may create new sections)
write_bodies = False            # choose whether to model bodies or not (experimental)
vortices_per_unit_length = 0.5  # resolution of vortex lattices
debug_geom = True
########################################################




# change directory and specify file name
loadpath = Path(filepath).absolute().parent
savepath = loadpath
path = Path(filepath).absolute()
os.chdir(path.parent)
filename = path.name
print('[VSP2AVL] Directory changed to "{}"'.format(path.parent))
print('[VSP2AVL] Now translating "{}"'.format(filename))

# airfoil_location = 'data/'
AVL_filename = re.split(r'\.',filename)[0] # sets AVL filename to same as DegenGeom file
component_delimiter = '# DegenGeom Type, Name, SurfNdx' # delimiter between component section of CSV


with open(filename) as f:
    DegenGeom = f.readlines()

components = []

# Look through CSV to find the starting line of each component. Save starting line index if it is a lifting surface, the name, and the ID of the component
for i, line in enumerate(DegenGeom):
    if component_delimiter in line:
        component = geom_data.geometry_component()
        component.begin_index = i

        next_line = re.split(r',\s*',DegenGeom[i+1])
        component.is_lifting_surface = next_line[0] == 'LIFTING_SURFACE'
        component.is_body = next_line[0] == 'BODY'

        component.name = next_line[1]
        component.num = next_line[2]
        component.ID = next_line[3]

        components.append(component)

# For each component, if it is a lifting surface, find the STICK_NODE section and save the beginning and end of that section
for i in range(len((components))):
    print(str(i) + ': ' + components[i].name + ' ' + str(components[i].num))
    # set ending index for component
    if i < len(components) - 1: #  if not last component in file
        components[i].end_index = components[i+1].begin_index # set ending index as beginning index of next component
    else: # if last component in file
        components[i].end_index = len(DegenGeom)

    components[i].find_SURFACE_NODE_data(DegenGeom)
    components[i].get_SURFACE_NODE_data(DegenGeom)
    components[i].get_body_profile()

    components[i].find_STICK_NODE_data(DegenGeom)
    if components[i].is_lifting_surface:
        components[i].get_STICK_NODE_data(DegenGeom)

    if component.is_lifting_surface and control_surfaces:
        components[i].get_control_surface_data(DegenGeom)
        components[i].interpret_control_surface(tolerance)

# calculate reference surface data
if Sref == 0 or Cref == 0 and refNum is not False:
    Sref = 0
    Cref = 0

    if components[refNum].stick_reoriented:
        for i in range(len(components[refNum].stick_le)-1):
            Sref += ((components[refNum].stick_chord[i]+components[refNum].stick_chord[i+1]) * components[refNum].stick_section_dist[i+1])/2
    else:
        for i in range(len(components[refNum].stick_le)-1):
            Sref += ((components[refNum].stick_chord[i]+components[refNum].stick_chord[i+1]) * components[refNum].stick_section_dist[i])/2

    Cref = Sref/np.sum(components[refNum].stick_section_dist)
    Bref = np.sum(components[refNum].stick_section_dist)

    Sref = Sref * 2
    Bref = Bref * 2


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
    if component.is_lifting_surface == True or (component.body_standard_body == False and write_bodies):
        AVL_file += component.create_lifting_surface(vortices_per_unit_length, loadpath, savepath)

    if component.is_body == True and component.body_standard_body == True and write_bodies:
        AVL_file += component.create_body()

# write AVL file array to new file
with open(AVL_filename + '.avl', 'w+') as f:
    for string in AVL_file:
        f.write(string)


if debug_geom:
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for component in components:
        if component.is_lifting_surface:
            le_x = []
            le_y = []
            le_z = []
            te_x = []
            te_y = []
            te_z = []
            for i in range(len(component.stick_le)):
                le_x.append(component.stick_le[i][0])
                le_y.append(component.stick_le[i][1])
                le_z.append(component.stick_le[i][2])
                te_x.append(le_x[i] + component.stick_chord[i]*np.cos(np.deg2rad(component.stick_Ainc[i])))
                te_y.append(le_y[i] + component.stick_chord[i]*np.sin(np.deg2rad(component.stick_Ainc[i]))*np.sin(component.stick_section_angle[i]))
                te_z.append(le_z[i] - component.stick_chord[i]*np.sin(np.deg2rad(component.stick_Ainc[i]))*np.cos(component.stick_section_angle[i]))

            if component.hingeline_data is not None:
                for j, name in enumerate(component.hingeline_name):
                    he_x = []
                    he_y = []
                    he_z = []
                    for i in range(len(component.stick_le)):
                        if component.hingeline_data[name]['is_here'][i]:
                            he_x.append(le_x[i] + component.stick_chord[i]*np.cos(np.deg2rad(component.stick_Ainc[i]))*component.hingeline_data[name]['x_c'][i])
                            he_y.append(le_y[i] + component.stick_chord[i]*np.sin(np.deg2rad(component.stick_Ainc[i]))*np.sin(component.stick_section_angle[i])*component.hingeline_data[name]['x_c'][i])
                            he_z.append(le_z[i] - component.stick_chord[i]*np.sin(np.deg2rad(component.stick_Ainc[i]))*np.cos(component.stick_section_angle[i])*component.hingeline_data[name]['x_c'][i])
                    ax.plot(he_x, he_y, he_z, color='orange')
                    # ax.scatter(he_x, he_y, he_z, color='orange')

            ax.plot(le_x, le_y, le_z, color='black')
            ax.plot(te_x, te_y, te_z, color='blue')
            # ax.scatter(le_x, le_y, le_z, color='black')
            # ax.scatter(te_x, te_y, te_z, color='blue')

            for i in range(len(le_x)):
                ax.plot([le_x[i], te_x[i]], [le_y[i], te_y[i]], [le_z[i], te_z[i]], color='green')

    ax.axis('equal')
    plt.show()
    

print(f"[VSP2AVL] AVL geometry file saved as \"{AVL_filename}.avl\"\n")