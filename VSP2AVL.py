import re
import os
import sys
import numpy as np
import geom_data
from pathlib import Path
import matplotlib.pyplot as plt
import configparser
import argparse


### IMPORT CONFIG DATA #######################################################
config = configparser.ConfigParser()

if Path('config.ini').is_file():
    config.read('config.ini')
else:
    config['AVL Parameters'] = {'use reference surface': 'False',
                                'reference surface Number': '1',
                                'sref': '0',
                                'cref': '0',
                                'bref': '0',
                                'xref': '0',
                                'yref': '0',
                                'zref': '0',
                                'mach number': '0'}
    
    config['VSP2AVL Settings'] = {'tolerance': '0.05',
                                  'model control surfaces': 'yes',
                                  'model bodies': 'no',
                                  'vortices per unit length': '0',
                                  'post processing geometry plot': 'no'}
    
    with open('config.ini', 'w') as configfile:
        config.write(configfile)
    print("Created 'config.ini' configuration file")

AVL_config = config['AVL Parameters']
VSP2AVL_config = config['VSP2AVL Settings']

useReferenceSurface = AVL_config.getboolean('use reference surface')
refNum = int(AVL_config['reference surface number'])
if not useReferenceSurface:
    refNum = False
Sref = float(AVL_config['sref'])                    # reference wing area
Cref = float(AVL_config['cref'])                    # reference chord length
Bref = float(AVL_config['bref'])                    # reference span
Xref = float(AVL_config['xref'])                    # center of gravity location
Yref = float(AVL_config['yref'])                    # center of gravity location
Zref = float(AVL_config['zref'])                    # center of gravity location
mach_number = float(AVL_config['mach number'])      # default mach number

tolerance = float(VSP2AVL_config['tolerance'])                                  # minimum allowed geometric distance between sections for hingeline section creation (make small)
control_surfaces = VSP2AVL_config.getboolean('model control surfaces')          # check for control surfaces in DegenGeom and add them to AVL file (may create new sections)
write_bodies = VSP2AVL_config.getboolean('model bodies')                        # choose whether to model bodies or not (experimental)
vortices_per_unit_length = float(VSP2AVL_config['vortices per unit length'])    # resolution of vortex lattices
debug_geom = VSP2AVL_config.getboolean('post processing geometry plot')
##############################################################################

parser = argparse.ArgumentParser(prog='VSP2AVL', description='Translate OpenVSP DegenGeom files into AVL files.')
parser.add_argument('filepath', type=str, nargs=1, help='DegenGeom CSV file path')
parser.add_argument('--debug', action='store_true', help='plot processed lifting surface data')
args = parser.parse_args()


filepath = os.path.join(args.filepath[0])
debug_geom += args.debug

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
        component = geom_data.geometry_component(DegenGeom,i)

        components.append(component)

if len(components) == 0:
    sys.exit('Error: Input file has no components. Ensure path to DegenGeom .csv is correct.')

# For each component, if it is a lifting surface, find the STICK_NODE section and save the beginning and end of that section
for i, component in enumerate(components):
    print(str(i) + ': ' + component.name + ' ' + str(component.num))
    # set ending index for component
    if i < len(components) - 1: #  if not last component in file
        components[i].end_index = components[i+1].begin_index # set ending index as beginning index of next component
    else: # if last component in file
        components[i].end_index = len(DegenGeom)

    components[i].get_SURFACE_NODE_data(DegenGeom)
    components[i].get_STICK_NODE_data(DegenGeom, tolerance)

    if component.is_lifting_surface and control_surfaces:
        components[i].get_control_surface_data(DegenGeom)
        components[i].interpret_control_surface(tolerance)

# calculate reference surface data
if refNum is not False and refNum is not None and len(components) != 0:
    Sref = 0
    Cref = 0

    if components[refNum].stick_reoriented:
        for i in range(len(components[refNum].stick_le[0])-1):
            Sref += ((components[refNum].stick_chord[0][i]+components[refNum].stick_chord[0][i+1]) * components[refNum].stick_section_dist[0][i+1])/2
    else:
        for i in range(len(components[refNum].stick_le[0])-1):
            Sref += ((components[refNum].stick_chord[0][i]+components[refNum].stick_chord[0][i+1]) * components[refNum].stick_section_dist[0][i])/2

    Cref = Sref/np.sum(components[refNum].stick_section_dist[0])
    Bref = np.sum(components[refNum].stick_section_dist[0])

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
    name_list = []
    hingeline_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
'#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
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
            for i in range(len(component.stick_le[0])):
                le_x.append(component.stick_le[0][i][0])
                le_y.append(component.stick_le[0][i][1])
                le_z.append(component.stick_le[0][i][2])
                te_x.append(le_x[i] + component.stick_chord[0][i]*np.cos(np.deg2rad(component.stick_Ainc[0][i])))
                te_y.append(le_y[i] + component.stick_chord[0][i]*np.sin(np.deg2rad(component.stick_Ainc[0][i]))*np.sin(component.stick_section_angle[0][i]))
                te_z.append(le_z[i] - component.stick_chord[0][i]*np.sin(np.deg2rad(component.stick_Ainc[0][i]))*np.cos(component.stick_section_angle[0][i]))

            if len(component.hingeline_data) != 0:
                for j, name in enumerate(component.hingeline_name):
                    if name not in name_list:
                        name_list.append(name)
                    color_num = np.mod(name_list.index(name),len(hingeline_colors))

                    he_x = []
                    he_y = []
                    he_z = []
                    for i in range(len(component.stick_le[0])):
                        if component.hingeline_data[name]['is_here'][i]:
                            he_x.append(le_x[i] + component.stick_chord[0][i]*np.cos(np.deg2rad(component.stick_Ainc[0][i]))*component.hingeline_data[name]['x_c'][i])
                            he_y.append(le_y[i] + component.stick_chord[0][i]*np.sin(np.deg2rad(component.stick_Ainc[0][i]))*np.sin(component.stick_section_angle[0][i])*component.hingeline_data[name]['x_c'][i])
                            he_z.append(le_z[i] - component.stick_chord[0][i]*np.sin(np.deg2rad(component.stick_Ainc[0][i]))*np.cos(component.stick_section_angle[0][i])*component.hingeline_data[name]['x_c'][i])
                    ax.plot(he_x, he_y, he_z, color=hingeline_colors[color_num])
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