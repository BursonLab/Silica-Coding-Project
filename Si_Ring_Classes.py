def xyz_to_list(file):
    """takes a file and converts it to a list"""
    text = [] #master list of all the points
    file_lines = file.readlines() #read all of the lines in the file
    file_lines = [x.strip() for x in file_lines] #strip off newline char
    file_lines = file_lines[2:] #remove header lines
    
    for line in file_lines:
        row = []  #create new line
        
        rowNum = int(line[0]) #find each chunk of data
        x = float(line[5:14])
        y = float(line[18:27])
        z = float(line[31:40])
        
        row.append(rowNum) #add to line
        row.append(x)
        row.append(y)
        row.append(z)
        text.append(row) #add line to master list
    return text


def get_distance(pt1, pt2):
    """ Finds the distance between two points. """
    x1, y1 = pt1
    x2, y2 = pt2
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** (1 / 2)


class Si():
    """ Contains the location of the Si atom, as well as each of the
        three rings surrounding it.  """

    """ Public Methods """

    def __init__(self, x, y, z):
        self._location = (x, y, z)
        self._rings = []
        self._d1 = 0
        self._d2 = 0
        self._d3 = 0

    def find_rings(self, ring_list):
        """ Finds the three rings bordering this Si atom, and stores
            them in self._rings. """
        self._findFirst(ring_list)
        if (len(self.get_rings() == 1)):
            self._findSecond()
        else:
            self._findThird(ring_list)

    def get_location(self):
        """ Returns the (x, y, z) of the atom. """
        return self._location

    def get_rings(self):
        """ Returns the list of rings bordering the atom. """
        return self._rings

    """ Private Methods """

    def _findFirst(self, ring_list):
        """ Finds the closest ring center to the atom. If there are
            equidistant centers, puts all into self._rings. """
        distance = 100000000000000000000
        answers = []
        for i in range(len(ring_list)):
            c1 = ring_list[i].get_location()
            c2 = self.get_location()
            if get_distance(c1, c2) < distance:
                answers.clear()
                answers.append(ring_list[i])
                distance = get_distance(c1, c2)
            if get_distance(c1, c2) == distance:
                answers.append(ring_list[i])
            for ring in answers:
                self._rings.append(ring)
                ring.set_atom(self)
            self._d1 = distance

    def _findSecond(self, ring_list):
        """ Finds the second closest ring center to the atom. If there
            are equidistant centers, puts all into self._rings. """
        distance = 100000000000000000000
        answers = []
        for i in range(len(ring_list)):
            c1 = ring_list[i].get_location()
            c2 = self.get_location()
            dist_2 = get_distance(c1, c2)
            if dist_2 < distance and dist_2 > self._d1:
                answers.clear()
                answers.append(ring_list[i])
                distance = dist_2
            if dist_2 == distance and dist_2 > self._d1:
                answers.append(ring_list[i])
            for ring in answers:
                self._rings.append(ring)
                ring.set_atom(self)
            self._d2 = distance

    def _findThird(self, ring_list):
        """ Finds the second closest ring center to the atom. """
        distance = 100000000000000000000
        answers = []
        for i in range(len(ring_list)):
            c1 = ring_list[i].get_location()
            c2 = self.get_location()
            dist_2 = get_distance(c1, c2)
            if dist_2 < distance and dist_2 > self._d2:
                answers.clear()
                answers.append(ring_list[i])
                distance = dist_2
            if dist_2 == distance and dist_2 > self._d2:
                answers.append(ring_list[i])
            for ring in answers:
                self._rings.append(ring)
                ring.set_atom(self)
            self._d3 = distance


class ring_center():
	""" Contains the location of the ring center, and the type of ring
        (number of members). """

    def __init__(self, x, y, z, ring_type):
        """ Constructor. """
        self._ring_type = ring_type
        self._location = (x, y, z)
        self._atoms = []

    def get_location(self):
    	""" Returns the location in (x, y, z) form. """
        	return self._location
    
    def get_type(self):
        """returns type of ring"""
        return self._ring_type

    def set_atom(self, atom);
        """ Puts an atom into self._atoms. """
        self._atoms.append(atom)
        
    def get_type(self):
        """returns type of ring"""
        return self._ring_type

    
def main():
    #get files
    si_file = open('Sample Si template.txt', encoding = 'utf-8')
    centers_file = open('Sample centers template.txt', encoding = 'utf-8')
    
    #convert to lists
    si_locs = xyz_to_list(si_file)
    centers = xyz_to_list(centers_file)

    for loc in si_locs:
        si = Si(loc[1], loc[2], loc[3])
        si.find_rings(centers)
        

main()
    
    
