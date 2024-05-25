### FUNCTIONS USED BY VSP2AVL
import numpy as np
import re

def get_STICK_NODE_data(component, lines, tolerance):
    '''
    Retrieve relevant data about component from STICK_NODE section of DegenGeom.
    '''
    stick_found = False
    component['stick_node_begin_index'] = []
    component['stick_node_end_index'] = []
    # iterate through related component lines to find the STICK_NODE section (used for AVL geometry)
    for j, line in enumerate(lines[component['begin_index']:component['end_index']]):
        if 'STICK_NODE' in line:
            component['stick_node_begin_index'].append(j)
            stick_found = True

        missing_end = len(component['stick_node_end_index']) < len(component['stick_node_begin_index'])

        # if STICK_NODE has been found for the component, find the next line that starts with #
        if stick_found and j > component['stick_node_begin_index'][-1] + 1 and line.startswith('#') and missing_end:
            component['stick_node_end_index'].append(j - 1)

    component['le'] = []
    component['te'] = []
    component['chord'] = []
    component['Ainc'] = []
    component['orig_num'] = []
    component['section_dist'] = []
    component['section_tolerance'] = []
    component['section_angle'] = []
    component['reoriented'] = False

    for m in range(len(component['stick_node_begin_index'])):
        # beginning index and ending index are set
        component['stick_node_begin_index'][m] += component['begin_index']
        component['stick_node_end_index'][m] += component['begin_index']

        # create new dictionary entries for AVL wing geometry
        component['le'].append([])
        component['te'].append([])
        component['chord'].append([])
        component['Ainc'].append([])
        component['orig_num'].append([])
        component['section_dist'].append([0])
        component['section_tolerance'].append([0])
        component['section_angle'].append([])

        # save data from STICK_NODE section
        for k,line in enumerate(lines[component['stick_node_begin_index'][m]+2:component['stick_node_end_index'][m]+1]):
            component['le'][m].append(np.float32(re.split(r',\s*',line)[:3])) # save LE location as tuple
            component['te'][m].append(np.float32(re.split(r',\s*',line)[3:6])) # save TE location as tuple
            
            # calculate chord length and angle of incidence from LE and TE location
            difference_vec = component['le'][m][-1] - component['te'][m][-1]
            chord_length = np.sqrt(np.sum(np.square(difference_vec))) # magnitude of difference vector
            component['chord'][m].append(chord_length)
            
            component['Ainc'][m].append(np.rad2deg(np.arcsin(np.dot(difference_vec, np.float32([0, 0, 1]))/chord_length)))
            # sin of dot product of difference vector and vector normal to XY plane ----> angle between difference vector and XY plane
            
            component['orig_num'][m].append(k)

            # calculate distance between sections, tolerance for each section, and angle between each section in YZ plane
            if k != 0:
                component['section_dist'][m].append(np.sqrt(np.sum(np.square(component['le'][m][k][1:3] - component['le'][m][k-1][1:3]))))
                component['section_tolerance'][m].append(component['section_dist'][m][k-1] * tolerance)
                component['section_angle'][m].append(np.arctan2(component['le'][m][k][2] - component['le'][m][k-1][2], component['le'][m][k][1] - component['le'][m][k-1][1]))
        component['section_dist'][m].append(0)
        component['section_tolerance'][m].append(0)
        component['section_angle'][m].insert(0,component['section_angle'][m][0])

        if 'standard_body' not in component or component['is_lifting_surface'] or component['standard_body']:
            # re-orient to spanwise low to high
            angle = np.arctan2(component['le'][m][1][2]-component['le'][m][0][2], component['le'][m][1][1]-component['le'][m][0][1])
            if quadrant(angle) == 2 or quadrant(angle) == 3:
                component['le'][m] = component['le'][m][::-1]
                component['te'][m] = component['te'][m][::-1]
                component['chord'][m] = component['chord'][m][::-1]
                component['Ainc'][m] = component['Ainc'][m][::-1]
                component['orig_num'][m] = component['orig_num'][m][::-1]
                component['section_dist'][m] = component['section_dist'][m][::-1]
                component['section_tolerance'][m] = component['section_tolerance'][m][::-1]
                component['section_angle'][m] = component['section_angle'][m][::-1]
                component['reoriented'] = True

        else:
            # re-orient to clockwise
            angle = np.arctan2(component['le'][m][1][2]-component['le'][m][0][2], component['le'][m][1][1]-component['le'][m][0][1])
            if quadrant(angle) == 1 or quadrant(angle) == 3:
                component['le'][m] = component['le'][m][::-1]
                component['te'][m] = component['te'][m][::-1]
                component['chord'][m] = component['chord'][m][::-1]
                component['Ainc'][m] = component['Ainc'][m][::-1]
                component['orig_num'][m] = component['orig_num'][m][::-1]
                component['section_dist'][m] = component['section_dist'][m][::-1]
                component['section_tolerance'][m] = component['section_tolerance'][m][::-1]
                component['section_angle'][m] = component['section_angle'][m][::-1]
                component['reoriented'] = True

    return component



def get_SURFACE_NODE_data(component, lines):
    surface_found = False
    # iterate through related component lines to find the SURFACE_NODE section (used for AVL geometry)
    for j, line in enumerate(lines[component['begin_index']:component['end_index']]):
        if 'SURFACE_NODE' in line:
            component['surface_node_begin_index'] = j
            component['num_secs'] = np.int16(re.split(r',\s*',line)[1])
            component['num_pts_per_sec'] = np.int16(re.split(r',\s*',line)[2])
            surface_found = True

        # if SURFACE_NODE has been found for the component, find the next line that starts with #
        if surface_found and j > component['surface_node_begin_index'] + 1 and line.startswith('SURFACE_FACE'):
            component['surface_node_end_index'] = j - 1
            break
            
    # beginning index and ending index are set
    component['surface_node_begin_index'] += component['begin_index']
    component['surface_node_end_index'] += component['begin_index']

    component['x_coords'] = []
    component['y_coords'] = []
    component['z_coords'] = []

    # save x, y, and z coordinates of surface nodes
    for k,line in enumerate(lines[component['surface_node_begin_index']+2:component['surface_node_end_index']+1]):
        component['x_coords'].append(np.float32(re.split(r',\s*',line)[0])) # save x coordinate
        component['y_coords'].append(np.float32(re.split(r',\s*',line)[1])) # save y coordinate
        component['z_coords'].append(np.float32(re.split(r',\s*',line)[2])) # save z coordinate

    x_max = np.max(component['x_coords'][:component['num_pts_per_sec']])
    x_min = np.min(component['x_coords'][:component['num_pts_per_sec']])

    if abs(x_max - x_min) > 1e-4:
        component['standard_body'] = False
    else:
        component['standard_body'] = True

    # get up and lo points for x-z profile of body
    component['up_coords'] = []
    component['lo_coords'] = []
    component['new_x_coords'] = []
    for s in range(component['num_secs']):
        up_point = True
        lo_point = True
        for p in range(component['num_pts_per_sec']):
            # print(s*components[i]['num_pts_per_sec'] + p)
            if up_point == True:
                up_point = component['z_coords'][s*component['num_pts_per_sec'] + p]
            else:
                if component['z_coords'][s*component['num_pts_per_sec'] + p] > up_point:
                    up_point = component['z_coords'][s*component['num_pts_per_sec'] + p]

            if lo_point == True:
                lo_point = component['z_coords'][s*component['num_pts_per_sec'] + p]
            else:
                if component['z_coords'][s*component['num_pts_per_sec'] + p] < lo_point:
                    lo_point = component['z_coords'][s*component['num_pts_per_sec'] + p]
        
        component['up_coords'].append(up_point)
        component['lo_coords'].append(lo_point)
        component['new_x_coords'].append(component['x_coords'][s*component['num_pts_per_sec']])
    
    # make common TE point for the body geometry if the last points aren't the same for the upper and lower surfaces
    if component['up_coords'][-1] != component['lo_coords'][-1]:
        mean = (component['up_coords'][-1] + component['lo_coords'][-1]) / 2
        component['up_coords'].append(mean)
        component['lo_coords'].append(mean)
        component['new_x_coords'].append(component['new_x_coords'][-1])

        # print(components[i]['up_coords'])
        # print(components[i]['lo_coords'])
        # print(components[i]['new_x_coords'])

    component['y-offset'] = np.average(component['y_coords'])

    reversed_up_coords = list(reversed(component['up_coords']))
    reversed_x_coords = list(reversed(component['new_x_coords']))

    component['body_y_points'] = reversed_up_coords + component['lo_coords'][1:]
    component['body_x_points'] = reversed_x_coords + component['new_x_coords'][1:]

    return component



def get_control_surface_data(component, lines):
    component['hingeline_index'] = []
    component['hingeline_name'] = []
    component['hingeline_start'] = []
    component['hingeline_end'] = []
    for i, line in enumerate(lines[component['begin_index']:component['end_index']]):
        if line.startswith('HINGELINE'):
            component['hingeline_index'].append(i+component['begin_index'])
            component['hingeline_name'].append(re.split(r',\s*',line)[1])
            component['hingeline_start'].append(np.array(np.float32(re.split(r',\s*',lines[i+component['begin_index']+2])[4:7])))
            component['hingeline_end'].append(np.array(np.float32(re.split(r',\s*',lines[i+component['begin_index']+2])[7:10])))

    # SORT HINGE START AND END IN SPANWISE DIRECTION FROM LOW TO HIGH
    if component['orig_num'][0][1] < component['orig_num'][0][0]: # (if component was reordered already)
        for i in range(len(component['hingeline_name'])):
            old_start = component['hingeline_start'][i]
            component['hingeline_start'][i] = component['hingeline_end'][i]
            component['hingeline_end'][i] = old_start

    return component



def interpret_control_surface(component, tolerance):
    if 'hingeline_name' in component and len(component['hingeline_name']) != 0:
        # print('has hingeline name')
        component['hingeline_data'] = {}

        for n, hingeline_name in enumerate(component['hingeline_name']):
            print("    ", hingeline_name)

            index_modifier = 0
            for orig_index in range(len(component['le'])-1):
                i = orig_index + index_modifier # index_modifier used if section is added due to hingeline start or end

                for iter in range(len(component['le'])-1):
                    component['section_dist'][iter+1] = np.sqrt(np.sum(np.square(component['le'][iter+1][1:3] - component['le'][iter][1:3])))
                    component['section_tolerance'][iter+1] = component['section_dist'][iter+1] * tolerance

                # print('current le: ', component['le'][i][1:3])
                # print('current hs: ', component['hingeline_start'][n][1:3])

                # project hinge difference vectors onto section difference vectors to check if hinge starts and ends are within tolerance of existing sections
                hinge_start_distance_vec = component['hingeline_start'][n][1:3] - component['le'][i][1:3]
                # print('current hsdv: ', hinge_start_distance_vec)
                hinge_start_distance = np.dot(hinge_start_distance_vec, component['le'][i+1][1:3] - component['le'][i][1:3]) / component['section_dist'][i+1]
                # print(hinge_start_distance)
                hinge_start_angle = np.arctan2(hinge_start_distance_vec[1], hinge_start_distance_vec[0])
                hinge_end_distance_vec = component['hingeline_end'][n][1:3] - component['le'][i][1:3]
                hinge_end_distance = np.dot(hinge_end_distance_vec, component['le'][i+1][1:3] - component['le'][i][1:3]) / component['section_dist'][i+1]
                hinge_end_angle = np.arctan2(hinge_end_distance_vec[1], hinge_end_distance_vec[0])
                #### implement angle tolerance check later #### TODO

                # check if hinge starts or ends at current section
                start_bounds = [hinge_start_distance - component['section_tolerance'][i+1], hinge_start_distance + component['section_tolerance'][i+1]]
                end_bounds = [hinge_end_distance - component['section_tolerance'][i+1], hinge_end_distance + component['section_tolerance'][i+1]]

                starts_here = start_bounds[0] < 0 and start_bounds[1] > 0
                ends_here = end_bounds[0] < 0 and end_bounds[1] > 0
                
                # print("starts here: ", starts_here)
                # print("ends here: ", ends_here)

                # check if control surface exists at current section
                is_here = hinge_start_distance < 0 and hinge_end_distance > 0
                # print("is here: ", is_here)

                # check if hinge starts or ends in the current section
                is_hingeline_start_between = hinge_start_distance > 0 and hinge_start_distance < (component['section_dist'][i+1] - component['section_tolerance'][i+1+1]) and not starts_here
                is_hingeline_end_between = hinge_end_distance > 0 and hinge_end_distance < (component['section_dist'][i+1] - component['section_tolerance'][i+1+1]) and not ends_here

                ### STARTS WITHIN?
                if is_hingeline_start_between:
                    # print('start between')
                    hingepoint = component['hingeline_start'][n]

                    # interpolate x, y, z coordinates, chord, and angle of incidence of section from sections on either side
                    x_dist_between_sects = component['le'][i+1][0]-component['le'][i][0] # get delta x
                    y_dist_between_sects = component['le'][i+1][1]-component['le'][i][1] # get delta y
                    z_dist_between_sects = component['le'][i+1][2]-component['le'][i][2] # get delta z
                    interped_x = component['le'][i][0] + (x_dist_between_sects/component['section_dist'][i+1])*hinge_start_distance
                    interped_y = component['le'][i][1] + (y_dist_between_sects/component['section_dist'][i+1])*hinge_start_distance
                    interped_z = component['le'][i][2] + (z_dist_between_sects/component['section_dist'][i+1])*hinge_start_distance
                    
                    hingeline_le = np.array([interped_x,interped_y,interped_z])
                    interped_chord = component['chord'][i] + ((component['chord'][i+1]-component['chord'][i])/component['section_dist'][i+1])*hinge_start_distance
                    interped_Ainc = component['Ainc'][i] + ((component['Ainc'][i+1]-component['Ainc'][i])/component['section_dist'][i+1])*hinge_start_distance
                    x_c_control_surface = abs(interped_x-hingepoint[0])/interped_chord

                    component['le'].insert(i+1,hingeline_le)
                    component['te'].insert(i+1,np.array([0,0,0]))
                    component['chord'].insert(i+1,interped_chord)
                    component['Ainc'].insert(i+1,interped_Ainc)
                    component['orig_num'].insert(i+1,(component['orig_num'][i+1]+component['orig_num'][i])/2)
                    component['section_dist'].insert(i+1+1,np.sqrt(np.sum(np.square(component['le'][i+1][1:3] - component['le'][i][1:3]))))
                    component['section_tolerance'].insert(i+1+1,component['section_dist'][i+1+1] * tolerance)

                    index_modifier += 1
                    i += 1

                for iter in range(len(component['le'])-1):
                    component['section_dist'][iter+1] = np.sqrt(np.sum(np.square(component['le'][iter+1][1:3] - component['le'][iter][1:3])))
                    component['section_tolerance'][iter+1] = component['section_dist'][iter+1] * tolerance

                ### ENDS WITHIN?
                if is_hingeline_end_between:
                    # print('end between')
                    hingepoint = component['hingeline_end'][n]

                    # interpolate x, y, z coordinates, chord, and angle of incidence of section from sections on either side
                    x_dist_between_sects = component['le'][i+1][0]-component['le'][i][0] # get delta x
                    y_dist_between_sects = component['le'][i+1][1]-component['le'][i][1] # get delta y
                    z_dist_between_sects = component['le'][i+1][2]-component['le'][i][2] # get delta z
                    interped_x = component['le'][i][0] + (x_dist_between_sects/component['section_dist'][i+1])*hinge_end_distance
                    interped_y = component['le'][i][1] + (y_dist_between_sects/component['section_dist'][i+1])*hinge_end_distance
                    interped_z = component['le'][i][2] + (z_dist_between_sects/component['section_dist'][i+1])*hinge_end_distance
                    
                    hingeline_le = np.array([interped_x,interped_y,interped_z])
                    interped_chord = component['chord'][i] + ((component['chord'][i+1]-component['chord'][i])/component['section_dist'][i+1])*hinge_end_distance
                    interped_Ainc = component['Ainc'][i] + ((component['Ainc'][i+1]-component['Ainc'][i])/component['section_dist'][i+1])*hinge_end_distance
                    x_c_control_surface = abs(interped_x-hingepoint[0])/interped_chord

                    component['le'].insert(i+1,hingeline_le)
                    component['te'].insert(i+1,np.array([0,0,0]))
                    component['chord'].insert(i+1,interped_chord)
                    component['Ainc'].insert(i+1,interped_Ainc)
                    component['orig_num'].insert(i+1,(component['orig_num'][i+1]+component['orig_num'][i])/2)
                    component['section_dist'].insert(i+1+1,np.sqrt(np.sum(np.square(component['le'][i+1][1:3] - component['le'][i][1:3]))))
                    component['section_tolerance'].insert(i+1+1,component['section_dist'][i+1+1] * tolerance)

                    index_modifier += 1

                for iter in range(len(component['le'])-1):
                    component['section_dist'][iter+1] = np.sqrt(np.sum(np.square(component['le'][iter+1][1:3] - component['le'][iter][1:3])))
                    component['section_tolerance'][iter+1] = component['section_dist'][iter+1] * tolerance

        # check to see if control surface should be defined for each section
        for n, hingeline_name in enumerate(component['hingeline_name']):

            component['hingeline_data'][hingeline_name] = {}
            component['hingeline_data'][hingeline_name]['is_here'] = [False] * len(component['le'][0])
            component['hingeline_data'][hingeline_name]['x_c'] = [0] * len(component['le'][0])


            for i in range(len(component['le'][0])-1):
                hinge_start_distance_vec = component['hingeline_start'][n][1:3] - component['le'][0][i][1:3]
                hinge_start_distance = np.dot(hinge_start_distance_vec, component['le'][0][i+1][1:3] - component['le'][0][i][1:3]) / component['section_dist'][0][i+1]
                hinge_end_distance_vec = component['hingeline_end'][n][1:3] - component['le'][0][i][1:3]
                hinge_end_distance = np.dot(hinge_end_distance_vec, component['le'][0][i+1][1:3] - component['le'][0][i][1:3]) / component['section_dist'][0][i+1]

                start_bounds = [hinge_start_distance - component['section_tolerance'][0][i+1], hinge_start_distance + component['section_tolerance'][0][i+1]]
                end_bounds = [hinge_end_distance - component['section_tolerance'][0][i+1], hinge_end_distance + component['section_tolerance'][0][i+1]]

                starts_here = start_bounds[0] < 0 and start_bounds[1] > 0
                ends_here = end_bounds[0] < 0 and end_bounds[1] > 0

                is_here = start_bounds[0] < 0 and hinge_end_distance - component['section_tolerance'][0][i+1] > 0

                # print(hinge_start_distance)
                # print(is_here)
                hinge_start_end_distance = np.sqrt(np.sum(np.square(component['hingeline_end'][n][1:3]-component['hingeline_start'][n][1:3])))
                hingeline_x_slope = (component['hingeline_end'][n][0]-component['hingeline_start'][n][0])/hinge_start_end_distance
                # print(hinge_start_end_distance)
                # print(hingeline_x_slope)
                # print(hinge_start_distance)
                interped_hingeline_x = component['hingeline_start'][n][0] + hingeline_x_slope * (-hinge_start_distance)
                x_c_control_surface = abs(interped_hingeline_x-component['le'][0][i][0])/component['chord'][0][i]

                component['hingeline_data'][hingeline_name]['x_c'][i] = x_c_control_surface

                # if control surface is here and does not need new starting section
                if is_here or starts_here or ends_here:
                    component['hingeline_data'][hingeline_name]['is_here'][i] = True

    return component


def create_lifting_surface(component, vortices_per_unit_length, loadpath, savepath):
    AVL_file = []
    surface_preamble = '''#--------------------------------------------------
SURFACE 
{} 
!Nchordwise  Cspace  Nspanwise  Sspace
12           1.0              

COMPONENT 
1

ANGLE
0.0

SCALE
1.0   1.0   1.0

TRANSLATE
0.0  0.0  0.0


'''.format(component['name']+' '+str(component['num']))
    AVL_file.append(surface_preamble)
    
    # create AVL section data
    for i in range(len(component['le'][0])):
        try:
            airfoil = '{}_{}_{}.dat'.format(component['name'], component['ID'], np.int16(np.round(component['orig_num'][0][i], 0)))

            # read airfoil data file
            with open(airfoil) as f:
                lines = f.readlines()

            # remove duplicate lines in airfoil data files
            lineset = set()
            j = 0
            while j < len(lines) - 1:
                if lines[j] in lineset:
                    lines.pop(j)
                    continue
                lineset.add(lines[j])
                j += 1

            # write airfoil data files changes
            with open(airfoil, 'w+') as f:
                f.writelines(lines)
            #print("Airfoil data file found: {}".format(airfoil))
        except FileNotFoundError:
            #print("No airfoil data found for section {} {} {}".format(component['name'], component['ID'], component['orig_num'][i]))
            airfoil = 'no_airfoil.dat'


        Xle, Yle, Zle = component['le'][0][i]
        chord = component['chord'][0][i]
        Ainc = component['Ainc'][0][i]
        Nspanwise = np.max([np.int16(vortices_per_unit_length*component['section_dist'][0][i+1]), 1])

        new_section = '''SECTION
#Xle    Yle    Zle     Chord   Ainc  Nspanwise  Sspace
{:.2f}    {:.2f}    {:.2f}     {:.2f}    {:.2f}    {}    0
'''.format(Xle, Yle, Zle, chord, Ainc, Nspanwise)

        #if isinstance(component['orig_num'][i], int):
        new_section += '''AFILE
{}
'''.format(airfoil)

        if 'hingeline_name' in component:
            for hingeline_name in component['hingeline_name']:
                if component['hingeline_data'][hingeline_name]['is_here'][i]:
                    new_section += '''CONTROL
{}     1.0   {:.2f}   0. 0. 0.   1   | name, gain,  Xhinge,  XYZhvec,  SgnDup
'''.format(hingeline_name, component['hingeline_data'][hingeline_name]['x_c'][i])

        new_section += '''
'''

        AVL_file.append(new_section) # append section data to AVL file array

    return AVL_file



def create_body(component):
    AVL_file = []
    body_geom_file_name = '{}_{}_{}.dat'.format(component['name'], component['ID'], component['num'])

    # create body geom data array
    geom_file = [body_geom_file_name + '\n']
    for i in range(len(component['body_x_points'])):
        geom_file.append('{} {}\n'.format(component['body_x_points'][i], component['body_y_points'][i]))

    body_entry = '''#=============================================
BODY
{}
28   1.0
#
TRANSLATE
0.0  {}  0.0
#
BFIL
{}

'''.format(component['name'] + ' ' + str(component['num']), component['y-offset'], body_geom_file_name)

    # write body data to new .dat file
    with open(body_geom_file_name, 'w+') as f:
        for string in geom_file:
            f.write(string)

    AVL_file.append(body_entry)

    return AVL_file



def quadrant(angle):
    if angle > np.pi:
        angle -= np.pi
    elif angle < -np.pi:
        angle += np.pi

    if angle < np.pi/2 and angle >= 0:
        quad = 1
    else:
        if angle <= np.pi and angle >= np.pi/2:
            quad = 2
        else:
            if angle <= -np.pi/2:
                quad = 3
            else:
                if angle < 0:
                    quad = 4

    return quad