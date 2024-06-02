### FUNCTIONS USED BY VSP2AVL
import numpy as np
import re

class geometry_component():

    def __init__(self):

        self.begin_index = None
        self.end_index = None
        
        self.is_lifting_surface = None
        self.is_body = None

        self.name = None
        self.num = None
        self.ID = None

        # STICK MODEL VARIABLES
        self.stick_node_begin_index = []
        self.stick_node_end_index = []
        self.stick_le = []
        self.stick_te = []
        self.stick_chord = []
        self.stick_Ainc = []
        self.stick_orig_num = []
        self.stick_section_dist = []
        self.stick_section_angle = []
        self.stick_reoriented = False

        # HINGELINE VARIABLES
        self.hingeline_index = []
        self.hingeline_name = []
        self.hingeline_start = []
        self.hingeline_end = []
        self.hingeline_data = {}

        # SURFACE DEGENGEOM DATA
        self.surface_node_begin_index = []
        self.surface_node_end_index = []
        self.surface_num_secs = None
        self.surface_num_pts_per_sec = None
        self.surface_coords = []

        # BODY MODEL VARIABLES
        self.body_standard_body = None
        self.body_up_coords = []
        self.body_lo_coords = []
        self.body_new_x_coords = []
        self.body_y_offset = None
        self.body_x_points = None
        self.body_y_points = None

    def find_STICK_NODE_data(self, lines):
        '''
        Find STICK_NODE starting and ending indicies.
        '''
        begin_index_found = False
        for j, line in enumerate(lines[self.begin_index:self.end_index]):
            if 'STICK_NODE' in line:
                self.stick_node_begin_index.append(j)
                begin_index_found = True

            if begin_index_found and j > self.stick_node_begin_index[-1] + 1 and line.startswith('#'):
                self.stick_node_end_index.append(j - 1)
                begin_index_found = False

        if len(self.stick_node_begin_index) == 1:
            self.stick_node_begin_index = self.stick_node_begin_index[0]
            self.stick_node_end_index = self.stick_node_end_index[0]

    def get_STICK_NODE_data(self, lines):
        '''
        Retrieve relevant data about component from STICK_NODE section of DegenGeom.
        '''
        self.stick_node_begin_index += self.begin_index
        self.stick_node_end_index += self.begin_index

        for k,line in enumerate(lines[self.stick_node_begin_index+2:self.stick_node_end_index+1]):
            self.stick_le.append(np.float32(re.split(r',\s*',line)[:3]))
            self.stick_te.append(np.float32(re.split(r',\s*',line)[3:6]))
            
            difference_vec = self.stick_le[-1] - self.stick_te[-1]
            chord_length = np.sqrt(np.sum(np.square(difference_vec)))
            
            self.stick_chord.append(chord_length)
            self.stick_Ainc.append(np.rad2deg(np.arcsin(np.dot(difference_vec, np.float32([0, 0, 1]))/chord_length)))
            self.stick_orig_num.append(k)

        self.stick_section_dist = np.zeros(len(self.stick_le)+1)
        self.stick_section_angle = np.zeros(len(self.stick_le)+1)

        if self.is_lifting_surface or self.body_standard_body:
            # reorder to spanwise low to high
            angle = np.arctan2(self.stick_le[1][2]-self.stick_le[0][2], self.stick_le[1][1]-self.stick_le[0][1])
            if self.quadrant(angle) == 2 or self.quadrant(angle) == 3:
                self.reorder_sections()

        else:
            # reorder to clockwise
            angle = np.arctan2(self.stick_le[1][2]-self.stick_le[0][2], self.stick_le[1][1]-self.stick_le[0][1])
            if self.quadrant(angle) == 1 or self.quadrant(angle) == 3:
                self.reorder_sections()

        self.recalc_sections()

    def reorder_sections(self):
        self.stick_le = self.stick_le[::-1]
        self.stick_te = self.stick_te[::-1]
        self.stick_chord = self.stick_chord[::-1]
        self.stick_Ainc = self.stick_Ainc[::-1]
        self.stick_orig_num = self.stick_orig_num[::-1]
        self.stick_section_dist = self.stick_section_dist[::-1]
        self.stick_section_angle = self.stick_section_angle[::-1]
        self.stick_reoriented = not self.stick_reoriented

    def recalc_sections(self):
        for iter in range(len(self.stick_le)-1):
            self.stick_section_dist[iter+1] = np.sqrt(np.sum(np.square(self.stick_le[iter+1][1:3] - self.stick_le[iter][1:3])))
            self.stick_section_angle[iter+1] = np.arctan2(self.stick_le[iter+1][2] - self.stick_le[iter][2], self.stick_le[iter+1][1] - self.stick_le[iter][1])

        
    def find_SURFACE_NODE_data(self, lines):
        begin_index_found = False
        for j, line in enumerate(lines[self.begin_index:self.end_index]):
            if 'SURFACE_NODE' in line:
                self.surface_node_begin_index.append(j)
                self.surface_num_secs = np.int16(re.split(r',\s*',line)[1])
                self.surface_num_pts_per_sec = np.int16(re.split(r',\s*',line)[2])
                begin_index_found = True

            if begin_index_found and j > self.surface_node_begin_index[-1] + 1 and line.startswith('SURFACE_FACE'):
                self.surface_node_end_index.append(j - 1)
                begin_index_found = False

        if len(self.surface_node_begin_index) == 1:
            self.surface_node_begin_index = self.surface_node_begin_index[0]
            self.surface_node_end_index = self.surface_node_end_index[0]

    def get_SURFACE_NODE_data(self, lines):
        '''
        Retrieve relevant data about component from SURFACE_NODE section of DegenGeom.
        '''
        self.surface_node_begin_index += self.begin_index
        self.surface_node_end_index += self.begin_index

        for k,line in enumerate(lines[self.surface_node_begin_index+2:self.surface_node_end_index+1]):
            self.surface_coords.append(np.float32(re.split(r',\s*',line)[:3]))

        plane_vec_1 = np.subtract(self.surface_coords[1], self.surface_coords[0])
        plane_vec_2 = np.subtract(self.surface_coords[2], self.surface_coords[0])
        plane_normal = np.cross(plane_vec_1, plane_vec_2)

        error = np.sum([np.dot(plane_normal, np.subtract(coord, self.surface_coords[0])) for coord in self.surface_coords[1:]])

        if abs(error) > 1e-4:
            self.body_standard_body = False
        else:
            self.body_standard_body = True

    def get_body_profile(self):
        for s in range(self.surface_num_secs):
            up_point = np.zeros(3)
            lo_point = np.zeros(3)
            for p in range(self.surface_num_pts_per_sec):
                if p == 0:
                    up_point = self.surface_coords[s*self.surface_num_pts_per_sec + p]
                elif self.surface_coords[s*self.surface_num_pts_per_sec + p][1] > up_point[1]:
                    up_point = self.surface_coords[s*self.surface_num_pts_per_sec + p]

                if p == 0:
                    lo_point = self.surface_coords[s*self.surface_num_pts_per_sec + p]
                elif self.surface_coords[s*self.surface_num_pts_per_sec + p][1] < lo_point[1]:
                    lo_point = self.surface_coords[s*self.surface_num_pts_per_sec + p]
            
            self.body_up_coords.append(up_point)
            self.body_lo_coords.append(lo_point)
            self.body_new_x_coords.append(self.surface_coords[s*self.surface_num_pts_per_sec][0])
        
        # make common TE point for the body geometry if the last points aren't the same for the upper and lower surfaces
        if self.body_up_coords[-1][1] != self.body_lo_coords[-1][1]:
            mean = np.add(self.body_up_coords[-1], self.body_lo_coords[-1]) / 2
            self.body_up_coords.append(mean)
            self.body_lo_coords.append(mean)
            self.body_new_x_coords.append(self.body_new_x_coords[-1])

        self.body_y_offset = np.average(self.surface_coords[:][1])

        reversed_up_coords = list(reversed(self.body_up_coords))
        reversed_x_coords = list(reversed(self.body_new_x_coords))

        self.body_y_points = reversed_up_coords + self.body_lo_coords[1:]
        self.body_x_points = reversed_x_coords + self.body_new_x_coords[1:]

        


    def get_control_surface_data(self, lines):
        for i, line in enumerate(lines[self.begin_index:self.end_index]):
            if line.startswith('HINGELINE'):
                self.hingeline_index.append(i+self.begin_index)
                self.hingeline_name.append(re.split(r',\s*',line)[1])
                self.hingeline_start.append(np.array(np.float32(re.split(r',\s*',lines[i+self.begin_index+2])[4:7])))
                self.hingeline_end.append(np.array(np.float32(re.split(r',\s*',lines[i+self.begin_index+2])[7:10])))

        if self.stick_reoriented:
            for i in range(len(self.hingeline_name)):
                [self.hingeline_start[i], self.hingeline_end[i]] = [self.hingeline_end[i], self.hingeline_start[i]]

    def interpret_control_surface(self, tolerance):
        if self.hingeline_name != None and len(self.hingeline_name) != 0:
            for n, hingeline_name in enumerate(self.hingeline_name):
                print("    ", hingeline_name)

                index_modifier = 0
                for orig_index in range(len(self.stick_le)-1):
                    i = orig_index + index_modifier # index_modifier used if section is added due to hingeline start or end

                    self.recalc_sections()

                    # project hinge difference vectors onto section difference vectors to check if hinge starts and ends are within tolerance of existing sections
                    hinge_start_distance_vec = self.hingeline_start[n][1:3] - self.stick_le[i][1:3]
                    hinge_start_distance = np.dot(hinge_start_distance_vec, np.subtract(self.stick_le[i+1][1:3], self.stick_le[i][1:3])) / self.stick_section_dist[i+1]
                    
                    hinge_start_angle = np.arctan2(hinge_start_distance_vec[1], hinge_start_distance_vec[0])
                    #### implement angle tolerance check later #### TODO

                    # check if hinge starts or ends at current section
                    start_bounds = [hinge_start_distance - tolerance, hinge_start_distance + tolerance]
                    starts_here = start_bounds[0] < 0 and start_bounds[1] > 0

                    # check if hinge starts in the current section
                    is_hingeline_start_between = hinge_start_distance >= tolerance and hinge_start_distance < (self.stick_section_dist[i+1] - tolerance) and not starts_here


                    ### STARTS WITHIN?
                    if is_hingeline_start_between:
                        # print('start between')
                        hingepoint = self.hingeline_start[n]

                        ## LEADING EDGE
                        # interpolate x, y, z coordinates from sections on either side
                        dist_between_sects = np.subtract(self.stick_le[i+1], self.stick_le[i])
                        interped_point = np.add(self.stick_le[i], (dist_between_sects/self.stick_section_dist[i+1]))*hinge_start_distance

                        ## TRAILING EDGE
                        # interpolate x, y, z coordinates, chord, and angle of incidence of section from sections on either side
                        te_x_dist_between_sects = self.stick_te[i+1][0]-self.stick_te[i][0] # get delta x
                        te_y_dist_between_sects = self.stick_te[i+1][1]-self.stick_te[i][1] # get delta y
                        te_z_dist_between_sects = self.stick_te[i+1][2]-self.stick_te[i][2] # get delta z
                        te_interped_x = self.stick_te[i][0] + (te_x_dist_between_sects/self.stick_section_dist[i+1])*hinge_start_distance
                        te_interped_y = self.stick_te[i][1] + (te_y_dist_between_sects/self.stick_section_dist[i+1])*hinge_start_distance
                        te_interped_z = self.stick_te[i][2] + (te_z_dist_between_sects/self.stick_section_dist[i+1])*hinge_start_distance
                        

                        hingeline_le = interped_point
                        hingeline_te = np.array([te_interped_x,te_interped_y,te_interped_z])
                        interped_chord = self.stick_chord[i] + ((self.stick_chord[i+1]-self.stick_chord[i])/self.stick_section_dist[i+1])*hinge_start_distance
                        interped_Ainc = self.stick_Ainc[i] + ((self.stick_Ainc[i+1]-self.stick_Ainc[i])/self.stick_section_dist[i+1])*hinge_start_distance
                        x_c_control_surface = abs(interped_point[0]-hingepoint[0])/interped_chord

                        self.stick_le.insert(i+1,hingeline_le)
                        self.stick_te.insert(i+1,hingeline_te)
                        self.stick_chord.insert(i+1,interped_chord)
                        self.stick_Ainc.insert(i+1,interped_Ainc)
                        self.stick_orig_num.insert(i+1,(self.stick_orig_num[i+1]+self.stick_orig_num[i])/2)
                        self.stick_section_dist = np.insert(self.stick_section_dist, i+1+1, np.sqrt(np.sum(np.square(np.subtract(self.stick_le[i+1][1:3], self.stick_le[i][1:3])))))
                        self.stick_section_angle = np.insert(self.stick_section_angle, i+1, (np.arctan2(self.stick_le[i+1][2] - self.stick_le[i][2], self.stick_le[i+1][1] - self.stick_le[i][1])))

                        index_modifier += 1
                        i += 1

                    self.recalc_sections()

                    # project hinge difference vectors onto section difference vectors to check if hinge starts and ends are within tolerance of existing sections
                    hinge_end_distance_vec = np.subtract(self.hingeline_end[n][1:3], self.stick_le[i][1:3])
                    hinge_end_distance = np.dot(hinge_end_distance_vec, np.subtract(self.stick_le[i+1][1:3], self.stick_le[i][1:3])) / self.stick_section_dist[i+1]
                    
                    hinge_end_angle = np.arctan2(hinge_end_distance_vec[1], hinge_end_distance_vec[0])
                    #### implement angle tolerance check later #### TODO
                    
                    # check if hinge starts or ends at current section
                    end_bounds = [hinge_end_distance - tolerance, hinge_end_distance + tolerance]
                    ends_here = end_bounds[0] < 0 and end_bounds[1] > 0

                    # check if hinge ends in the current section
                    is_hingeline_end_between = hinge_end_distance >= tolerance and hinge_end_distance < (self.stick_section_dist[i+1] - tolerance) and not ends_here

                    ### ENDS WITHIN?
                    if is_hingeline_end_between:
                        # print('end between')
                        hingepoint = self.hingeline_end[n]

                        ## LEADING EDGE
                        # interpolate x, y, z coordinates, chord, and angle of incidence of section from sections on either side
                        x_dist_between_sects = self.stick_le[i+1][0]-self.stick_le[i][0] # get delta x
                        y_dist_between_sects = self.stick_le[i+1][1]-self.stick_le[i][1] # get delta y
                        z_dist_between_sects = self.stick_le[i+1][2]-self.stick_le[i][2] # get delta z
                        interped_x = self.stick_le[i][0] + (x_dist_between_sects/self.stick_section_dist[i+1])*hinge_end_distance
                        interped_y = self.stick_le[i][1] + (y_dist_between_sects/self.stick_section_dist[i+1])*hinge_end_distance
                        interped_z = self.stick_le[i][2] + (z_dist_between_sects/self.stick_section_dist[i+1])*hinge_end_distance

                        ## TRAILING EDGE
                        # interpolate x, y, z coordinates, chord, and angle of incidence of section from sections on either side
                        te_x_dist_between_sects = self.stick_te[i+1][0]-self.stick_te[i][0] # get delta x
                        te_y_dist_between_sects = self.stick_te[i+1][1]-self.stick_te[i][1] # get delta y
                        te_z_dist_between_sects = self.stick_te[i+1][2]-self.stick_te[i][2] # get delta z
                        te_interped_x = self.stick_te[i][0] + (te_x_dist_between_sects/self.stick_section_dist[i+1])*hinge_end_distance
                        te_interped_y = self.stick_te[i][1] + (te_y_dist_between_sects/self.stick_section_dist[i+1])*hinge_end_distance
                        te_interped_z = self.stick_te[i][2] + (te_z_dist_between_sects/self.stick_section_dist[i+1])*hinge_end_distance

                        
                        hingeline_le = np.array([interped_x,interped_y,interped_z])
                        hingeline_te = np.array([te_interped_x,te_interped_y,te_interped_z])
                        interped_chord = self.stick_chord[i] + ((self.stick_chord[i+1]-self.stick_chord[i])/self.stick_section_dist[i+1])*hinge_end_distance
                        interped_Ainc = self.stick_Ainc[i] + ((self.stick_Ainc[i+1]-self.stick_Ainc[i])/self.stick_section_dist[i+1])*hinge_end_distance
                        x_c_control_surface = abs(interped_x-hingepoint[0])/interped_chord

                        self.stick_le.insert(i+1,hingeline_le)
                        self.stick_te.insert(i+1,hingeline_te)
                        self.stick_chord.insert(i+1,interped_chord)
                        self.stick_Ainc.insert(i+1,interped_Ainc)
                        self.stick_orig_num.insert(i+1,(self.stick_orig_num[i+1]+self.stick_orig_num[i])/2)
                        self.stick_section_dist = np.insert(self.stick_section_dist, i+1+1,np.sqrt(np.sum(np.square(self.stick_le[i+1][1:3] - self.stick_le[i][1:3]))))
                        self.stick_section_angle = np.insert(self.stick_section_angle, i+1, (np.arctan2(self.stick_le[i+1][2] - self.stick_le[i][2], self.stick_le[i+1][1] - self.stick_le[i][1])))

                        index_modifier += 1

                    self.recalc_sections()

            # check to see if control surface should be defined for each section
            for n, hingeline_name in enumerate(self.hingeline_name):

                self.hingeline_data[hingeline_name] = {}
                self.hingeline_data[hingeline_name]['is_here'] = [False] * len(self.stick_le)
                self.hingeline_data[hingeline_name]['x_c'] = [0] * len(self.stick_le)


                for i in range(len(self.stick_le)-1):
                    hinge_start_distance_vec = self.hingeline_start[n][1:3] - self.stick_le[i][1:3]
                    hinge_start_distance = np.dot(hinge_start_distance_vec, self.stick_le[i+1][1:3] - self.stick_le[i][1:3]) / self.stick_section_dist[i+1]
                    hinge_end_distance_vec = self.hingeline_end[n][1:3] - self.stick_le[i][1:3]
                    hinge_end_distance = np.dot(hinge_end_distance_vec, self.stick_le[i+1][1:3] - self.stick_le[i][1:3]) / self.stick_section_dist[i+1]

                    start_bounds = [hinge_start_distance - tolerance, hinge_start_distance + tolerance]
                    end_bounds = [hinge_end_distance - tolerance, hinge_end_distance + tolerance]

                    starts_here = start_bounds[0] < 0 and start_bounds[1] > 0
                    ends_here = end_bounds[0] < 0 and end_bounds[1] > 0

                    is_here = start_bounds[0] < 0 and hinge_end_distance > 0


                    hinge_start_end_distance = np.sqrt(np.sum(np.square(self.hingeline_end[n][1:3]-self.hingeline_start[n][1:3])))
                    hingeline_x_slope = (self.hingeline_end[n][0]-self.hingeline_start[n][0])/hinge_start_end_distance

                    interped_hingeline_x = self.hingeline_start[n][0] + hingeline_x_slope * (-hinge_start_distance)
                    x_c_control_surface = abs(interped_hingeline_x-self.stick_le[i][0])/self.stick_chord[i]

                    self.hingeline_data[hingeline_name]['x_c'][i] = x_c_control_surface

                    # if control surface is here and does not need new starting section
                    if is_here or starts_here or ends_here:
                        self.hingeline_data[hingeline_name]['is_here'][i] = True

        


    def create_lifting_surface(self, vortices_per_unit_length, loadpath, savepath):
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


'''.format(self.name+' '+str(self.num))
        AVL_file.append(surface_preamble)
        
        # create AVL section data
        for i in range(len(self.stick_le[0])):
            try:
                airfoil = '{}_{}_{}.dat'.format(self.name, self.ID, np.int16(np.round(self.stick_orig_num[i], 0)))

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

            except FileNotFoundError:
                airfoil = 'no_airfoil.dat'


            Xle, Yle, Zle = self.stick_le[i]
            chord = self.stick_chord[i]
            Ainc = self.stick_Ainc[i]
            Nspanwise = np.max([np.int16(np.ceil(vortices_per_unit_length*self.stick_section_dist[i+1])), 1])

            new_section = '''SECTION
#Xle    Yle    Zle     Chord   Ainc  Nspanwise  Sspace
{:.2f}    {:.2f}    {:.2f}     {:.2f}    {:.2f}    {}    0
'''.format(Xle, Yle, Zle, chord, Ainc, Nspanwise)

            #if isinstance(self.orig_num'][i], int):
            new_section += '''AFILE
{}
'''.format(airfoil)

            if self.hingeline_name != None:
                for hingeline_name in self.hingeline_name:
                    if self.hingeline_data[hingeline_name]['is_here'][i]:
                        new_section += '''CONTROL
{}     1.0   {:.2f}   0. 0. 0.   1   | name, gain,  Xhinge,  XYZhvec,  SgnDup
'''.format(hingeline_name, self.hingeline_data[hingeline_name]['x_c'][i])

            new_section += '''
'''

            AVL_file.append(new_section) # append section data to AVL file array

        return AVL_file


    def create_body(self):
        AVL_file = []
        body_geom_file_name = '{}_{}_{}.dat'.format(self.name, self.ID, self.num)

        # create body geom data array
        geom_file = [body_geom_file_name + '\n']
        for i in range(len(self.body_x_points)):
            geom_file.append('{} {}\n'.format(self.body_x_points[i], self.body_y_points[i]))

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

'''.format(self.name + ' ' + str(self.num), self.body_y_offset, body_geom_file_name)

        # write body data to new .dat file
        with open(body_geom_file_name, 'w+') as f:
            for string in geom_file:
                f.write(string)

        AVL_file.append(body_entry)

        return AVL_file


    @staticmethod
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