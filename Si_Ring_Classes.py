from sklearn.neighbors import NearestNeighbors


def get_distance(pt1, pt2):
    """ Finds the distance between two points. """
    x1 = pt1[0]
    y1 = pt1[1]
    x2 = pt2[0]
    y2 = pt2[1]
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** (1 / 2)


class Si():
    """ Contains the location of the Si atom, as well as each of the
        three rings surrounding it.  """

    """ Public Methods """

    def __init__(self, x, y, z):
        self._location = [x, y, z]
        self._rings = []
        self._d1 = 0
        self._d2 = 0
        self._d3 = 0

    def _findClosestThree(self, ring_list, x_max, y_max, edge_buffer):
        if self.is_edge(x_max, y_max, edge_buffer):
            return
        
        ring_pos = []
        
        
        for ring in ring_list:
            ring_pos.append(ring.get_location())
            
        nearest = NearestNeighbors(n_neighbors=3, algorithm='ball_tree').fit(ring_pos)
        dist, ind = nearest.kneighbors([self.get_location()])
        for i in range(len(ind[0])):
            self._rings.append(ring_list[ind[0][i]])
        
        

    
    def find_rings(self, ring_list, x_max, y_max, edge_buffer):
        """ Finds the three rings bordering this Si atom, and stores
            them in self._rings.
        self._findFirst(ring_list, x_max, y_max, edge_buffer)
        print("1st ring found!")
        if (len(self.get_rings()) == 1):
            self._findSecond(ring_list, x_max, y_max, edge_buffer)
        if (len(self.get_rings()) == 2):
            self._findThird(ring_list, x_max, y_max, edge_buffer)
        print("size: ", len(self._rings))"""
        self._findClosestThree(ring_list, x_max, y_max, edge_buffer)

    def get_location(self):
        """ Returns the (x, y, z) of the atom. """
        return self._location

    def get_rings(self):
        """ Returns the list of rings bordering the atom. """
        return self._rings

    def is_edge(self, max_x, max_y, edge_buffer):
        """ Determines if this Si atom is on the edge of the image
            returns true if so, false otherwise. """
        x = self.get_location()[0]
        y = self.get_location()[1]
        d = edge_buffer
        return x < d or x > max_x - d or y < d or y > max_y - d

    """ Private Methods """

    def _findFirst(self, ring_list, x_max, y_max, edge_buffer):
        """ Finds the closest ring center to the atom. If there are
            equidistant centers, puts all into self._rings. """
        # Excludes any Si atoms that are included as an edge case
        if self.is_edge(x_max, y_max, edge_buffer):
            return
        
        # Sets an arbitrary number as the first distance. This number
        # is used because it will be bigger than any distance
        # calculated.
        distance = 100000000000000000000
        answers = []
        for i in range(len(ring_list)):
            c1 = ring_list[i].get_location()
            c2 = self.get_location()
            
            # Checks if the calculate distance is less than the current
            # smallest distance. If so, resets the answer list and adds
            # the newest ring.
            if get_distance(c1, c2) < distance:
                answers = []
                answers.append(ring_list[i])
                distance = get_distance(c1, c2)
            if get_distance(c1, c2) == distance:
                answers.append(ring_list[i])
            for ring in answers:
                self._rings.append(ring)
                # print(len(self._rings))
                ring.set_atom(self)
            self._d1 = distance

    def _findSecond(self, ring_list, x_max, y_max, edge_buffer):
        """ Finds the second closest ring center to the atom. If there
            are equidistant centers, puts all into self._rings. """
        if self.is_edge(x_max, y_max, edge_buffer):
            return
        distance = 100000000000000000000
        answers = []
        for i in range(len(ring_list)):
            c1 = ring_list[i].get_location()
            c2 = self.get_location()
            dist_2 = get_distance(c1, c2)
            if dist_2 < distance and dist_2 > self._d1:
                answers = []
                answers.append(ring_list[i])
                distance = dist_2
            if dist_2 == distance and dist_2 > self._d1:
                answers.append(ring_list[i])
            for ring in answers:
                self._rings.append(ring)
                ring.set_atom(self)
            self._d2 = distance

    def _findThird(self, ring_list, x_max, y_max, edge_buffer):
        """ Finds the second closest ring center to the atom. """
        if self.is_edge(x_max, y_max, edge_buffer):
            return
        distance = 100000000000000000000
        answers = []
        for i in range(len(ring_list)):
            c1 = ring_list[i].get_location()
            c2 = self.get_location()
            dist_2 = get_distance(c1, c2)
            if dist_2 < distance and dist_2 > self._d2:
                answers = []
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
        (number of members). Objects should always be constructed using the nm 
        coordinates.  Nanometers is the assumed unit."""

    def __init__(self, ring_type, x, y, z, unit):
        """ Constructor. """
        self._ring_type = ring_type
        if unit == "nm":
            self._nm_location = [x, y, z]
            self._pixel_location = [0, 0, 0]
        else:
            self._nm_location = [0, 0, 0]
            self._pixel_location = [x, y, z]
        self._atoms = []

    def get_nm_location(self):
        """ Returns the location in (x, y, z) form. Units are nm. """
        return self._nm_location
    
    def get_pix_location(self):
        """ Returns the location in (x, y, z) form. Units are Pixels"""
        return self._location
    
    def find_nm_location(self, nm_dim, im_width, im_height):
        """ Finds the coordinates in nm when the pixel coordinates are 
        known. """
        scale = (nm_dim[0]/im_width)
        for i in range(3):
            self._nm_location[i] = scale * self._pixels_location[i]

    def find_pix_location(self, nm_dim, im_width, im_height):
        """ Finds the coordinates in pixels when the nm coordinates are 
        known. """
        scale = (im_width/nm_dim[0])
        for i in range(3):
            self._pixel_location[i] = scale * self._location[i]

    def change_location(self, x, y, z, unit, nm_dim, im_width, im_height):
        """ Changes the coordinates of the center, and finds the coordinates in
        the other unit. """
        if unit == "nm":
            self._nm_location = [x, y, z]
            self.find_pix_location(nm_dim, im_width, im_height)
        else:
            self._pixel_location = [x, y, z]
            self.find_nm_location(nm_dim, im_width, im_height)

    def get_type(self):
        """returns type of ring"""
        return self._ring_type

    def set_atom(self, atom):
        """ Puts an atom into self._atoms. """
        self._atoms.append(atom)

    def get_atoms(self):
        """ Returns the atom list """
        return self._atoms

    def remove(self, index):
        """ Removes an atom from the atom list BY INDEX """
        del self._atoms[index]


def main():
    return
main()
