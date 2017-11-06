def xyz_to_list(file):
    """takes a file and converts it to a list"""
    text = []  # master list of all the points
    file_lines = file.readlines()  # read all of the lines in the file
    file_lines = [x.strip() for x in file_lines]  # strip off newline char
    file_lines = file_lines[2:]  # remove header lines

    for line in file_lines:
        row = []  # create new line

        rowNum = int(line[0])  # find each chunk of data
        x = float(line[5:14])
        y = float(line[18:27])
        z = float(line[31:40])

        row.append(rowNum)  # add to line
        row.append(x)
        row.append(y)
        row.append(z)
        text.append(row)  # add line to master list
    return text


def order(lst):
    """ Returns a new list with the original's data, sorted smallest to
        largest. Specifically for use with lists of length 3. """
    ordered = []
    smallest = lst[0]
    while len(lst) != 0:
        for i in range(len(lst)):
            if lst[0] < smallest:
                smallest = lst[0]
        ordered.append(smallest)
        lst.remove(smallest)
    return ordered


def find_type(atom):
    """ Determines the type of an Si atom's triplet. Returns that type
        in smallest-largest order. """
    rings = atom.get_rings()
    t1 = rings[0].get_type()
    t2 = rings[1].get_type()
    t3 = rings[2].get_type()
    ordered = order([t1, t2, t3])
    return (ordered[0], ordered[1], ordered[2])


def is_present(lst, target):
    """ Determines of the target is in the lst. If so, returns true. If
        not, returns false. """
    for item in lst:
        if target == item:
            return True
    return False


def get_stats(si_list):
    """ Determines the number of each type of ring triplet [(5, 5, 6),
        (5, 6, 7), etc] and returns a list of tuples and ints containing the
        triplet type and the number found: [(5, 6, 7), 10, (5, 5, 5), 15,
        etc]. """
    types = []
    for atom in si_list:
        typ = find_type(atom)
        if is_present(types, typ):
            types[types.index(typ) + 1] += 1
        else:
            types.append(typ)
            types.append(1)
    return types


def get_distance(pt1, pt2):
    """ Finds the distance between two points. """
    x1, y1, z1 = pt1
    x2, y2, z2 = pt2
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) ** (1 / 2)


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

    def find_rings(self, ring_list, x_max, y_max, edge_buffer):
        """ Finds the three rings bordering this Si atom, and stores
            them in self._rings. """
        self._findFirst(ring_list, x_max, y_max, edge_buffer)
        if (len(self.get_rings()) == 1):
            self._findSecond(ring_list, x_max, y_max, edge_buffer)
        if (len(self.get_rings()) == 2):
            self._findThird(ring_list, x_max, y_max, edge_buffer)

    def get_location(self):
        """ Returns the (x, y, z) of the atom. """
        return self._location

    def get_rings(self):
        """ Returns the list of rings bordering the atom. """
        return self._rings

    def is_edge(self, max_x, max_y, edge_buffer):
        """ Determines if this Si atom is on the edge of the image
            returns true if so, false otherwise. """

        # Uses arbitrary value to determine edge cases (10 currently)
        x, y, _ = self.get_location()
        d = edge_buffer
        return x < d or x > max_x - d or y < d or y > max_y - d

    """ Private Methods """

    def _findFirst(self, ring_list, x_max, y_max, edge_buffer):
        """ Finds the closest ring center to the atom. If there are
            equidistant centers, puts all into self._rings. """
        if self.is_edge(x_max, y_max, edge_buffer):
            return
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
                answers.clear()
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

    def __init__(self, x, y, z):  # removed , ring_type from arguments
        """ Constructor. """
#        self._ring_type = ring_type
        self._location = (x, y, z)
        self._atoms = []

    def get_location(self):
        """ Returns the location in (x, y, z) form. """
        return self._location

#    def get_type(self):
#        """returns type of ring"""
#        return self._ring_type

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
#    # get files
#    si_file = open('Sample Si template.txt', encoding='utf-8')
#    centers_file = open('Sample centers template.txt', encoding='utf-8')
#
#    # convert to lists
#    si_locs = xyz_to_list(si_file)
#    centers = xyz_to_list(centers_file)
#
#    for loc in si_locs:
#        si = Si(loc[1], loc[2], loc[3])
#        si.find_rings(centers)
#
#    for si in si_locs:
#        print(si.get_location(), end=" ")
#        for ring in si.get_rings():
#            print(ring.get_location(), end=" ")
#        print()


main()
