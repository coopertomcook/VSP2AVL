### FUNCTIONS USED BY VSP2AVL
import numpy as np
import re

class geometry_component():

    def __init__(self, file, begin_index):
        self.end_index = None

        self.begin_index = begin_index
        next_line = re.split(r',\s*', file[begin_index+1])

        self.is_lifting_surface = next_line[0] == 'LIFTING_SURFACE'
        self.is_body = next_line[0] == 'BODY'

        self.name = next_line[1]
        self.num = next_line[2]
        self.ID = next_line[3]

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
        self.surface_node_begin_index = None
        self.surface_node_end_index = None
        self.surface_num_secs = None
        self.surface_num_pts_per_sec = None
        self.surface_x_coords = []
        self.surface_y_coords = []
        self.surface_z_coords = []

        # BODY MODEL VARIABLES
        self.body_standard_body = None
        self.body_up_coords = []
        self.body_lo_coords = []
        self.body_new_x_coords = []
        self.body_y_offset = None
        self.body_x_points = None
        self.body_y_points = None


    def get_STICK_NODE_data(self, lines, tolerance):
        '''
        Retrieve relevant data about component from STICK_NODE section of DegenGeom.
        '''
        stick_found = False
        # iterate through related component lines to find the STICK_NODE section (used for AVL geometry)
        for j, line in enumerate(lines[self.begin_index:self.end_index]):
            if 'STICK_NODE' in line:
                self.stick_node_begin_index.append(j)
                stick_found = True

            missing_end = len(self.stick_node_end_index) < len(self.stick_node_begin_index)

            # if STICK_NODE has been found for the component, find the next line that starts with #
            if stick_found and j > self.stick_node_begin_index[-1] + 1 and line.startswith('#') and missing_end:
                self.stick_node_end_index.append(j - 1)

        for m in range(len(self.stick_node_begin_index)):
            # beginning index and ending index are set
            self.stick_node_begin_index[m] += self.begin_index
            self.stick_node_end_index[m] += self.begin_index

            # create new dictionary entries for AVL wing geometry
            self.stick_le.append([])
            self.stick_te.append([])
            self.stick_chord.append([])
            self.stick_Ainc.append([])
            self.stick_orig_num.append([])
            self.stick_section_dist.append([0])
            self.stick_section_angle.append([0])

            # save data from STICK_NODE section
            for k,line in enumerate(lines[self.stick_node_begin_index[m]+2:self.stick_node_end_index[m]+1]):
                self.stick_le[m].append(np.float32(re.split(r',\s*',line)[:3])) # save LE location as tuple
                self.stick_te[m].append(np.float32(re.split(r',\s*',line)[3:6])) # save TE location as tuple
                
                # calculate chord length and angle of incidence from LE and TE location
                difference_vec = self.stick_le[m][-1] - self.stick_te[m][-1]
                chord_length = np.sqrt(np.sum(np.square(difference_vec))) # magnitude of difference vector
                self.stick_chord[m].append(chord_length)
                
                self.stick_Ainc[m].append(np.rad2deg(np.arcsin(np.dot(difference_vec, np.float32([0, 0, 1]))/chord_length)))
                # sin of dot product of difference vector and vector normal to XY plane ----> angle between difference vector and XY plane
                
                self.stick_orig_num[m].append(k)

                # calculate distance between sections, tolerance for each section, and angle between each section in YZ plane
                if k != 0:
                    self.stick_section_dist[m].append(np.sqrt(np.sum(np.square(self.stick_le[m][k][1:3] - self.stick_le[m][k-1][1:3]))))
                    self.stick_section_angle[m].append(np.arctan2(self.stick_le[m][k][2] - self.stick_le[m][k-1][2], self.stick_le[m][k][1] - self.stick_le[m][k-1][1]))
            self.stick_section_dist[m].append(0)
            self.stick_section_angle[m].append(0)

            if self.body_standard_body == None or self.is_lifting_surface or self.body_standard_body:
                # re-orient to spanwise low to high
                angle = np.arctan2(self.stick_le[m][1][2]-self.stick_le[m][0][2], self.stick_le[m][1][1]-self.stick_le[m][0][1])
                if self.quadrant(angle) == 2 or self.quadrant(angle) == 3:
                    self.stick_le[m] = self.stick_le[m][::-1]
                    self.stick_te[m] = self.stick_te[m][::-1]
                    self.stick_chord[m] = self.stick_chord[m][::-1]
                    self.stick_Ainc[m] = self.stick_Ainc[m][::-1]
                    self.stick_orig_num[m] = self.stick_orig_num[m][::-1]
                    self.stick_section_dist[m] = self.stick_section_dist[m][::-1]
                    self.stick_section_angle[m] = self.stick_section_angle[m][::-1]
                    self.stick_reoriented = True

            else:
                # re-orient to clockwise
                angle = np.arctan2(self.stick_le[m][1][2]-self.stick_le[m][0][2], self.stick_le[m][1][1]-self.stick_le[m][0][1])
                if self.quadrant(angle) == 1 or self.quadrant(angle) == 3:
                    self.stick_le[m] = self.stick_le[m][::-1]
                    self.stick_te[m] = self.stick_te[m][::-1]
                    self.stick_chord[m] = self.stick_chord[m][::-1]
                    self.stick_Ainc[m] = self.stick_Ainc[m][::-1]
                    self.stick_orig_num[m] = self.stick_orig_num[m][::-1]
                    self.stick_section_dist[m] = self.stick_section_dist[m][::-1]
                    self.stick_section_angle[m] = self.stick_section_angle[m][::-1]
                    self.stick_reoriented = True

        self.recalc_sections()

        


    def get_SURFACE_NODE_data(self, lines):
        surface_found = False
        # iterate through related component lines to find the SURFACE_NODE section (used for AVL geometry)
        for j, line in enumerate(lines[self.begin_index:self.end_index]):
            if 'SURFACE_NODE' in line:
                self.surface_node_begin_index = j
                self.surface_num_secs = np.int16(re.split(r',\s*',line)[1])
                self.surface_num_pts_per_sec = np.int16(re.split(r',\s*',line)[2])
                surface_found = True

            # if SURFACE_NODE has been found for the component, find the next line that starts with #
            if surface_found and j > self.surface_node_begin_index + 1 and line.startswith('SURFACE_FACE'):
                self.surface_node_end_index = j - 1
                break
                
        # beginning index and ending index are set
        self.surface_node_begin_index += self.begin_index
        self.surface_node_end_index += self.begin_index

        # save x, y, and z coordinates of surface nodes
        for k,line in enumerate(lines[self.surface_node_begin_index+2:self.surface_node_end_index+1]):
            self.surface_x_coords.append(np.float32(re.split(r',\s*',line)[0])) # save x coordinate
            self.surface_y_coords.append(np.float32(re.split(r',\s*',line)[1])) # save y coordinate
            self.surface_z_coords.append(np.float32(re.split(r',\s*',line)[2])) # save z coordinate

        x_max = np.max(self.surface_x_coords[:self.surface_num_pts_per_sec])
        x_min = np.min(self.surface_x_coords[:self.surface_num_pts_per_sec])

        if abs(x_max - x_min) > 1e-4:
            self.body_standard_body = False
        else:
            self.body_standard_body = True

        # get up and lo points for x-z profile of body
        for s in range(self.surface_num_secs):
            up_point = True
            lo_point = True
            for p in range(self.surface_num_pts_per_sec):
                # print(s*components[i]['num_pts_per_sec'] + p)
                if up_point == True:
                    up_point = self.surface_z_coords[s*self.surface_num_pts_per_sec + p]
                else:
                    if self.surface_z_coords[s*self.surface_num_pts_per_sec + p] > up_point:
                        up_point = self.surface_z_coords[s*self.surface_num_pts_per_sec + p]

                if lo_point == True:
                    lo_point = self.surface_z_coords[s*self.surface_num_pts_per_sec + p]
                else:
                    if self.surface_z_coords[s*self.surface_num_pts_per_sec + p] < lo_point:
                        lo_point = self.surface_z_coords[s*self.surface_num_pts_per_sec + p]
            
            self.body_up_coords.append(up_point)
            self.body_lo_coords.append(lo_point)
            self.body_new_x_coords.append(self.surface_x_coords[s*self.surface_num_pts_per_sec])
        
        # make common TE point for the body geometry if the last points aren't the same for the upper and lower surfaces
        if self.body_up_coords[-1] != self.body_lo_coords[-1]:
            mean = (self.body_up_coords[-1] + self.body_lo_coords[-1]) / 2
            self.body_up_coords.append(mean)
            self.body_lo_coords.append(mean)
            self.body_new_x_coords.append(self.body_new_x_coords[-1])

        self.body_y_offset = np.average(self.surface_y_coords)

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

        # SORT HINGE START AND END IN SPANWISE DIRECTION FROM LOW TO HIGH
        if self.stick_reoriented: # (if component was reordered already)
            for i in range(len(self.hingeline_name)):
                old_start = self.hingeline_start[i]
                self.hingeline_start[i] = self.hingeline_end[i]
                self.hingeline_end[i] = old_start

    def recalc_sections(self):
        for iter in range(len(self.stick_le[0])-1):
            self.stick_section_dist[0][iter+1] = np.sqrt(np.sum(np.square(self.stick_le[0][iter+1][1:3] - self.stick_le[0][iter][1:3])))
            self.stick_section_angle[0][iter+1] = np.arctan2(self.stick_le[0][iter+1][2] - self.stick_le[0][iter][2], self.stick_le[0][iter+1][1] - self.stick_le[0][iter][1])


    def interpret_control_surface(self, tolerance):
        if self.hingeline_name != None and len(self.hingeline_name) != 0:
            for n, hingeline_name in enumerate(self.hingeline_name):
                print("    ", hingeline_name)

                index_modifier = 0
                for orig_index in range(len(self.stick_le[0])-1):
                    i = orig_index + index_modifier # index_modifier used if section is added due to hingeline start or end

                    self.recalc_sections()

                    # project hinge difference vectors onto section difference vectors to check if hinge starts and ends are within tolerance of existing sections
                    hinge_start_distance_vec = self.hingeline_start[n][1:3] - self.stick_le[0][i][1:3]
                    hinge_start_distance = np.dot(hinge_start_distance_vec, self.stick_le[0][i+1][1:3] - self.stick_le[0][i][1:3]) / self.stick_section_dist[0][i+1]
                    
                    hinge_start_angle = np.arctan2(hinge_start_distance_vec[1], hinge_start_distance_vec[0])
                    #### implement angle tolerance check later #### TODO

                    # check if hinge starts or ends at current section
                    start_bounds = [hinge_start_distance - tolerance, hinge_start_distance + tolerance]
                    starts_here = start_bounds[0] < 0 and start_bounds[1] > 0

                    # check if control surface exists at current section
                    is_here = hinge_start_distance < 0 and hinge_end_distance > 0

                    # check if hinge starts in the current section
                    is_hingeline_start_between = hinge_start_distance >= tolerance and hinge_start_distance < (self.stick_section_dist[0][i+1] - tolerance) and not starts_here
                    

                    # print('current le: ', self.stick_le[0][i][1:3])
                    # print('current hs: ', self.hingeline_start[n][1:3])
                    # print('start_vec:  ', hinge_start_distance_vec)
                    # print('current he: ', self.hingeline_end[n][1:3])
                    # print('end_vec:    ', hinge_end_distance_vec)

                    # print(f"here  = {is_here}")
                    # print(f"start = {is_hingeline_start_between}")
                    # print(f"end   = {is_hingeline_end_between}")

                    ### STARTS WITHIN?
                    if is_hingeline_start_between:
                        # print('start between')
                        hingepoint = self.hingeline_start[n]

                        ## LEADING EDGE
                        # interpolate x, y, z coordinates, chord, and angle of incidence of section from sections on either side
                        x_dist_between_sects = self.stick_le[0][i+1][0]-self.stick_le[0][i][0] # get delta x
                        y_dist_between_sects = self.stick_le[0][i+1][1]-self.stick_le[0][i][1] # get delta y
                        z_dist_between_sects = self.stick_le[0][i+1][2]-self.stick_le[0][i][2] # get delta z
                        interped_x = self.stick_le[0][i][0] + (x_dist_between_sects/self.stick_section_dist[0][i+1])*hinge_start_distance
                        interped_y = self.stick_le[0][i][1] + (y_dist_between_sects/self.stick_section_dist[0][i+1])*hinge_start_distance
                        interped_z = self.stick_le[0][i][2] + (z_dist_between_sects/self.stick_section_dist[0][i+1])*hinge_start_distance

                        ## TRAILING EDGE
                        # interpolate x, y, z coordinates, chord, and angle of incidence of section from sections on either side
                        te_x_dist_between_sects = self.stick_te[0][i+1][0]-self.stick_te[0][i][0] # get delta x
                        te_y_dist_between_sects = self.stick_te[0][i+1][1]-self.stick_te[0][i][1] # get delta y
                        te_z_dist_between_sects = self.stick_te[0][i+1][2]-self.stick_te[0][i][2] # get delta z
                        te_interped_x = self.stick_te[0][i][0] + (te_x_dist_between_sects/self.stick_section_dist[0][i+1])*hinge_start_distance
                        te_interped_y = self.stick_te[0][i][1] + (te_y_dist_between_sects/self.stick_section_dist[0][i+1])*hinge_start_distance
                        te_interped_z = self.stick_te[0][i][2] + (te_z_dist_between_sects/self.stick_section_dist[0][i+1])*hinge_start_distance
                        

                        hingeline_le = np.array([interped_x,interped_y,interped_z])
                        hingeline_te = np.array([te_interped_x,te_interped_y,te_interped_z])
                        interped_chord = self.stick_chord[0][i] + ((self.stick_chord[0][i+1]-self.stick_chord[0][i])/self.stick_section_dist[0][i+1])*hinge_start_distance
                        interped_Ainc = self.stick_Ainc[0][i] + ((self.stick_Ainc[0][i+1]-self.stick_Ainc[0][i])/self.stick_section_dist[0][i+1])*hinge_start_distance
                        x_c_control_surface = abs(interped_x-hingepoint[0])/interped_chord

                        self.stick_le[0].insert(i+1,hingeline_le)
                        self.stick_te[0].insert(i+1,hingeline_te)
                        self.stick_chord[0].insert(i+1,interped_chord)
                        self.stick_Ainc[0].insert(i+1,interped_Ainc)
                        self.stick_orig_num[0].insert(i+1,(self.stick_orig_num[0][i+1]+self.stick_orig_num[0][i])/2)
                        self.stick_section_dist[0].insert(i+1+1,np.sqrt(np.sum(np.square(self.stick_le[0][i+1][1:3] - self.stick_le[0][i][1:3]))))
                        self.stick_section_angle[0].insert(i+1, (np.arctan2(self.stick_le[0][i+1][2] - self.stick_le[0][i][2], self.stick_le[0][i+1][1] - self.stick_le[0][i][1])))

                        index_modifier += 1
                        i += 1

                    self.recalc_sections()

                    # project hinge difference vectors onto section difference vectors to check if hinge starts and ends are within tolerance of existing sections
                    hinge_end_distance_vec = self.hingeline_end[n][1:3] - self.stick_le[0][i][1:3]
                    hinge_end_distance = np.dot(hinge_end_distance_vec, self.stick_le[0][i+1][1:3] - self.stick_le[0][i][1:3]) / self.stick_section_dist[0][i+1]
                    
                    hinge_end_angle = np.arctan2(hinge_end_distance_vec[1], hinge_end_distance_vec[0])
                    #### implement angle tolerance check later #### TODO
                    
                    # check if hinge starts or ends at current section
                    end_bounds = [hinge_end_distance - tolerance, hinge_end_distance + tolerance]
                    ends_here = end_bounds[0] < 0 and end_bounds[1] > 0
                    
                    # check if control surface exists at current section
                    is_here = hinge_start_distance < 0 and hinge_end_distance > 0
                    
                    # check if hinge ends in the current section
                    is_hingeline_end_between = hinge_end_distance >= tolerance and hinge_end_distance < (self.stick_section_dist[0][i+1] - tolerance) and not ends_here

                    ### ENDS WITHIN?
                    if is_hingeline_end_between:
                        # print('end between')
                        hingepoint = self.hingeline_end[n]

                        ## LEADING EDGE
                        # interpolate x, y, z coordinates, chord, and angle of incidence of section from sections on either side
                        x_dist_between_sects = self.stick_le[0][i+1][0]-self.stick_le[0][i][0] # get delta x
                        y_dist_between_sects = self.stick_le[0][i+1][1]-self.stick_le[0][i][1] # get delta y
                        z_dist_between_sects = self.stick_le[0][i+1][2]-self.stick_le[0][i][2] # get delta z
                        interped_x = self.stick_le[0][i][0] + (x_dist_between_sects/self.stick_section_dist[0][i+1])*hinge_end_distance
                        interped_y = self.stick_le[0][i][1] + (y_dist_between_sects/self.stick_section_dist[0][i+1])*hinge_end_distance
                        interped_z = self.stick_le[0][i][2] + (z_dist_between_sects/self.stick_section_dist[0][i+1])*hinge_end_distance

                        ## TRAILING EDGE
                        # interpolate x, y, z coordinates, chord, and angle of incidence of section from sections on either side
                        te_x_dist_between_sects = self.stick_te[0][i+1][0]-self.stick_te[0][i][0] # get delta x
                        te_y_dist_between_sects = self.stick_te[0][i+1][1]-self.stick_te[0][i][1] # get delta y
                        te_z_dist_between_sects = self.stick_te[0][i+1][2]-self.stick_te[0][i][2] # get delta z
                        te_interped_x = self.stick_te[0][i][0] + (te_x_dist_between_sects/self.stick_section_dist[0][i+1])*hinge_end_distance
                        te_interped_y = self.stick_te[0][i][1] + (te_y_dist_between_sects/self.stick_section_dist[0][i+1])*hinge_end_distance
                        te_interped_z = self.stick_te[0][i][2] + (te_z_dist_between_sects/self.stick_section_dist[0][i+1])*hinge_end_distance

                        
                        hingeline_le = np.array([interped_x,interped_y,interped_z])
                        hingeline_te = np.array([te_interped_x,te_interped_y,te_interped_z])
                        interped_chord = self.stick_chord[0][i] + ((self.stick_chord[0][i+1]-self.stick_chord[0][i])/self.stick_section_dist[0][i+1])*hinge_end_distance
                        interped_Ainc = self.stick_Ainc[0][i] + ((self.stick_Ainc[0][i+1]-self.stick_Ainc[0][i])/self.stick_section_dist[0][i+1])*hinge_end_distance
                        x_c_control_surface = abs(interped_x-hingepoint[0])/interped_chord

                        self.stick_le[0].insert(i+1,hingeline_le)
                        self.stick_te[0].insert(i+1,hingeline_te)
                        self.stick_chord[0].insert(i+1,interped_chord)
                        self.stick_Ainc[0].insert(i+1,interped_Ainc)
                        self.stick_orig_num[0].insert(i+1,(self.stick_orig_num[0][i+1]+self.stick_orig_num[0][i])/2)
                        self.stick_section_dist[0].insert(i+1+1,np.sqrt(np.sum(np.square(self.stick_le[0][i+1][1:3] - self.stick_le[0][i][1:3]))))
                        self.stick_section_angle[0].insert(i+1, (np.arctan2(self.stick_le[0][i+1][2] - self.stick_le[0][i][2], self.stick_le[0][i+1][1] - self.stick_le[0][i][1])))

                        index_modifier += 1

                    self.recalc_sections()

            # check to see if control surface should be defined for each section
            for n, hingeline_name in enumerate(self.hingeline_name):

                self.hingeline_data[hingeline_name] = {}
                self.hingeline_data[hingeline_name]['is_here'] = [False] * len(self.stick_le[0])
                self.hingeline_data[hingeline_name]['x_c'] = [0] * len(self.stick_le[0])


                for i in range(len(self.stick_le[0])-1):
                    hinge_start_distance_vec = self.hingeline_start[n][1:3] - self.stick_le[0][i][1:3]
                    hinge_start_distance = np.dot(hinge_start_distance_vec, self.stick_le[0][i+1][1:3] - self.stick_le[0][i][1:3]) / self.stick_section_dist[0][i+1]
                    hinge_end_distance_vec = self.hingeline_end[n][1:3] - self.stick_le[0][i][1:3]
                    hinge_end_distance = np.dot(hinge_end_distance_vec, self.stick_le[0][i+1][1:3] - self.stick_le[0][i][1:3]) / self.stick_section_dist[0][i+1]

                    start_bounds = [hinge_start_distance - tolerance, hinge_start_distance + tolerance]
                    end_bounds = [hinge_end_distance - tolerance, hinge_end_distance + tolerance]

                    starts_here = start_bounds[0] < 0 and start_bounds[1] > 0
                    ends_here = end_bounds[0] < 0 and end_bounds[1] > 0

                    is_here = start_bounds[0] < 0 and hinge_end_distance > 0

                    # print(hinge_start_distance)
                    # print(is_here)
                    hinge_start_end_distance = np.sqrt(np.sum(np.square(self.hingeline_end[n][1:3]-self.hingeline_start[n][1:3])))
                    hingeline_x_slope = (self.hingeline_end[n][0]-self.hingeline_start[n][0])/hinge_start_end_distance
                    # print(hinge_start_end_distance)
                    # print(hingeline_x_slope)
                    # print(hinge_start_distance)
                    interped_hingeline_x = self.hingeline_start[n][0] + hingeline_x_slope * (-hinge_start_distance)
                    x_c_control_surface = abs(interped_hingeline_x-self.stick_le[0][i][0])/self.stick_chord[0][i]

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
                airfoil = '{}_{}_{}.dat'.format(self.name, self.ID, np.int16(np.round(self.stick_orig_num[0][i], 0)))

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


            Xle, Yle, Zle = self.stick_le[0][i]
            chord = self.stick_chord[0][i]
            Ainc = self.stick_Ainc[0][i]
            Nspanwise = np.max([np.int16(np.ceil(vortices_per_unit_length*self.stick_section_dist[0][i+1])), 1])

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
            angle -= 2*np.pi
        elif angle < -np.pi:
            angle += 2*np.pi

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