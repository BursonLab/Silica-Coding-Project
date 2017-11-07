import math
import numpy
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import Si_Ring_Classes


def distance(position1, position2):
    """finds the distance between two atoms"""
    return math.sqrt(math.pow(position1[0] - position2[0], 2) +
                     math.pow(position1[1] - position2[1], 2) +
                     math.pow(position1[2] - position2[2], 2))


def dists(positions, dist):
    """finds if a triplet could have an Si atom between them"""

    # if there were not enough close to make a triplet, return none
    if len(positions) < 3:
        return[""]
    # if there is a triplet and they are close enough to have a Si,
    # return the triplet, else return blank
    if len(positions) == 3:
        if distance(positions[1], positions[2]) <= dist:
            return positions
        else:
            return[""]
    numbers = []

#    if len(positions) == 5:
#        print(1)
    # if there are more then 2 close enough to have a Si between them, findthe
    #  one that could not given the other two
    for i in range(len(positions)):
        numbers.append(0)
    for i in range(1, len(positions) - 1):
        for j in range(1, len(positions) - i):
            # if two positions are not close enough, add a counter to both.
            # If they are close enough, remove a counter from both
            if distance(positions[i], positions[i + j]) > dist:
                numbers[i] += 1
                numbers[i + j] += 1
            else:
                numbers[i] -= 1
                numbers[i + j] -= 1

    # removetheonewiththemostcounters
    del positions[numbers.index(max(numbers))]

    # if these still are not close enough to have a triplet between them,
    # return none. If they are close enough, return the new triplet
    if distance(positions[1], positions[2]) <= dist:
        return positions
    else:
        return[""]


def find_four(opositions, far):
    """finds four membered rings and returns a list of lists of their
       locaitons"""
    rings = [[]]
    remov = []
    # for each oxygen
    for i in range(len(opositions)):
        rings.append([""])
        rings[i] = [opositions[i]]
        # for each oxygen with an x position higher than the current
        for j in range(1, len(opositions) - i):
            # if th exposition is less than the possible distance between two
            #  oxygenatoms(variableinclusionradius)
            if abs(opositions[i][0] - opositions[i + j][0]) <= far:
                # if the distance between the two oxygens is less than the
                #  characteristic distance(variable inclusion radius)
                if distance(opositions[i], opositions[i + j]) <= far:
                    rings[i].append(opositions[i + j])
        rem = 0
        if len(rings[i]) < 4:
            rem = 1
        elif len(rings[i]) > 4:
            while len(rings[i]) != 4:
                distances = []
                for k in range(len(rings[i])):
                    tot_len = 0
                    for l in range(1, len(rings[i]) - k):
                        tot_len += distance(rings[i][k], rings[i][k + l])
                    distances.append(tot_len)
                del rings[i][distances.index(max(distances))]
        if len(rings[i]) == 4:
            distances = []
            for n in range(len(rings[i]) - 1):
                for m in range(1, len(rings[i]) - n):
                    distances.append(distance(rings[i][n], rings[i][n + m]))
            for n in range(2):
                del distances[distances.index(max(distances))]
            for n in range(4):
                for m in range(1, len(distances) - n):
                    if abs(distances[n] - distances[n + m]) > .03:
                        rem = 1
        if rem == 1:
            remov.insert(0, i)

    for n in range(len(remov)):
        del rings[remov[n]]

    return rings


def triarea(p1, p2, p3):
    """finds the area of triangle"""
    a = distance(p1, p2)
    b = distance(p2, p3)
    c = distance(p1, p3)
    s = (a + b + c) / 2
    return math.sqrt(s * (s - a) * (s - b) * (s - c))


def ringarea(corners):
    n = len(corners)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return float(area)


def rem4(rings, si):
    """finds if the silicon atom is within a 4 membered ring"""
    for i in range(len(rings)):
        triangles = 0
        distances = []
        locations = []
        for n in range(len(rings[i]) - 1):
            for m in range(1, len(rings[i]) - n):
                distances.append(distance(rings[i][n], rings[i][n + m]))
                locations.append([n, n + m])
        locations.append(len(rings[i]))
        for n in range(2):
            del locations[distances.index(max(distances))]
            del distances[distances.index(max(distances))]
        for n in range(len(locations)):
            triangles += triarea(rings[i][locations[n][0]],
                                 rings[i][locations[n][1]], si)
        if ringarea(rings[i]) == triangles:
            return"n"
    return"y"


def si_finder(opositions):
    """finds the position of a Si given a triplet of oxygen"""

    # characteristic distance
    dist = 1.6 * math.pow(10, - 1)

    # sets up the translation to happen around a basepoint(the first point in
    #  the positions)
    trans = [[0, 0, 0], [opositions[1][0] - opositions[0][0],
                         opositions[1][1] - opositions[0][1],
                         opositions[1][2] - opositions[0][2]],
             [opositions[2][0] - opositions[0][0],
              opositions[2][1] - opositions[0][1],
              opositions[2][2] - opositions[0][2]]]

    # finds vector perpendicular to the plane of the three points
    v = numpy.matrix([numpy.linalg.det([[trans[1][1], trans[2][1]],
                                        [trans[1][2], trans[2][2]]]),
                      numpy.linalg.det([[trans[1][0], trans[2][0]],
                                        [trans[1][2], trans[2][2]]]),
                      numpy.linalg.det([[trans[1][0], trans[2][0]],
                                        [trans[1][1], trans[2][1]]])])

    # sets up first rotation matrix about the x axis
    theta = math.atan2(v.item(1), v.item(2))
    xmatr = numpy.matrix([[1, 0, 0], [0, math.cos(theta), - math.sin(theta)],
                          [0, math.sin(theta), math.cos(theta)]])
    trans1 = numpy.matrix(trans)
    rot1 = numpy.matrix.dot(trans1, xmatr)
    v1 = numpy.matrix.dot(v, xmatr)

    # second rotation matrix about the y axis
    rho = math.atan2(v1.item(0), v1.item(2))
    ymatr = numpy.matrix([[math.cos(rho), 0, math.sin(rho)], [0, 1, 0],
                          [-math.sin(rho), 0, math.cos(rho)]])
    rot2 = numpy.matrix.dot(rot1, ymatr)

    # should be in the xy plane now. Have to rotate such that two points
    #  are on the x axis
    alph = math.atan2(rot2.item(4), rot2.item(3))
    bet = math.atan2(rot2.item(7), rot2.item(6))
    r1 = math.sqrt(math.pow(rot2.item(3), 2) + math.pow(rot2.item(4), 2))
    r2 = math.sqrt(math.pow(rot2.item(6), 2) + math.pow(rot2.item(7), 2))
    x = r1 / 2
    y = r2 * (1 - math.cos(bet - alph)) / (2.0 * math.sin(bet - alph))
    z = math.sqrt(abs(math.pow(dist, 2) - math.pow(x, 2) - math.pow(y, 2)))
    si_pos = numpy.matrix([x, y, z])

    # rotate back to originial position
    init = math.atan2(si_pos.item(1), si_pos.item(0))
    r = math.sqrt(math.pow(si_pos.item(0), 2) + math.pow(si_pos.item(1), 2))
    x = r * math.cos(init + alph)
    y = r * math.sin(init + alph)
    si_pos = numpy.matrix([x, y, z])

    # undo second rotation matrix
    iymatr = numpy.linalg.inv(ymatr)
    si_pos = numpy.matrix.dot(si_pos, iymatr)

    # undo first rotation matrix
    ixmatr = numpy.linalg.inv(xmatr)
    si_pos = numpy.matrix.dot(si_pos, ixmatr)

    # translate back so there is no point at the origin
    si_pos = [si_pos.item(0) + opositions[0][0],
              si_pos.item(1) + opositions[0][1],
              si_pos.item(2) + opositions[0][2]]

    return si_pos


def o_locator(opositions):
    """locates all possiblee triplets"""

    # assumed oxygens are ordered by increasing x values
    # used to collect all the found oxygens close enough to have a single Si
    #  between them
    found = [[""]]
    # for each oxygen
    for i in range(len(opositions)):
        found[i] = [opositions[i]]
        # for each oxygen with an x position higher than the current
        for j in range(1, len(opositions) - i):
            # if the x position is less than the possible distance between two
            #  oxygenatoms(variableinclusionradius)
            if abs(opositions[i][0] - opositions[i + j][0]) <= \
                    3.45 * math.pow(10, - 1):
                # if the distance between the two oxygens is less than the
                #  characteristic distance(variable inclusion radius)
                if distance(opositions[i], opositions[i + j]) <= \
                        3.45 * math.pow(10, - 1):
                    found[i].append(opositions[i + j])
        found.append([""])

    # removes last appended empty list
    del found[len(found) - 1]

    # remove all those too far apart using dist function (variable inclusion
    #  radius)
    for n in range(len(found)):
        found[n] = dists(found[n], .345)

    # createanarrayforpositionstoremove
    remov = []
    # for all atoms with found oxygens
    for n in range(len(found)):
        # add empties to a list for removal
        if found[n] == [""]:
            remov.insert(0, n)

    # remove those in the remove list
    for m in range(len(remov)):
        del found[remov[m]]

    # return the list of those oxygen that have a possible Si between them
    return found


def locate_si(positions, dist):
    # assumes presorted positions by x position
    doubles = []

    # finds all within the given radius and adds those doubles to the list
    for i in range(len(positions)):
        for j in range(1, len(positions) - i):
            if distance(positions[i], positions[i + j]) <= dist:
                doubles.append([positions[i], positions[i + j]])

    return doubles


def find_o(positions, dist):

    opositions = []

    for i in range(len(positions)):
        # center at origin
        pos2 = [positions[i][1][0] - positions[i][0][0], positions[i][1][1] -
                positions[i][0][1], positions[i][1][2] - positions[i][0][2]]

        # rotate until both points are in the xy plane
        theta = numpy.arctan2(pos2[1], pos2[0])
        phi = numpy.arctan2(pos2[2], pos2[0])
        newx = math.sqrt(math.pow(pos2[0], 2) + math.pow(pos2[2], 2))
        newy = newx * math.tan(theta)

        # find in si position (midpoint between origin and pos 2 in the x - y
        #  plane with x making up the difference)
        x = newx / 2
        y = newy / 2

        if math.pow(dist, 2) - math.pow(x, 2) - math.pow(y, 2) > 0:
            z = math.sqrt(math.pow(dist, 2) - math.pow(x, 2) - math.pow(y, 2))
        else:
            z = 0

        # current angle above x - y plane
        r = math.sqrt(math.pow(x, 2) + math.pow(y, 2) + math.pow(z, 2))
        alph = math.asin(z / r)

        # when rotated back, it will rotate to angle phi + alph
        opos = [r * math.cos(theta) * math.cos(alph + phi),
                r * math.sin(theta) * math.cos(alph + phi),
                r * math.sin(alph + phi)]

        # append to the list
        opositions.append([opos[0] + positions[i][0][0], opos[1] +
                           positions[i][0][1], opos[2] + positions[i][0][2]])

    return opositions

def read_xyz_line(line):
    row = []
    val = ""
    
    for digit in line:
        digit = str(digit)
        if (digit == " "):
            val = float(val)
            row.append(val)
            val = ""
        else:
            val += digit       
    return row

def xyz_to_list(file):
    text = []
    with open(file) as f:
        file_lines = f.readline()
        while(file_lines):
            row = read_xyz_line(file_lines)
            text.append(row)
            file_lines = f.readline()
#            print(row)
    return text

def xyz_to_objects(file):
    center_xyz_list = xyz_to_list(file)
    center_obj_list = []
    for c in center_xyz_list:
        center = Si_Ring_Classes.ring_center(c[0], c[1], c[2], 0)
        center_obj_list.append(center)
    return center_obj_list

#stat functions
def order(lst):
    """ Returns a new list with the original's data, sorted smallest to
        largest. """
    ordered = []
    while len(lst) != 0:
        smallest = lst[0]
        for i in range(len(lst)):
            if lst[i] < smallest:
                smallest = lst[i]
        ordered.append(smallest)
        lst.remove(smallest)
    return ordered


def find_type(atom):
    """ Determines the type of an Si atom's triplet. Returns that type
        in smallest-largest order. """
    rings = atom.get_rings()
    print(rings)
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




def main():

    # input center positions
    cpfile = "CPSampleXYZ.txt" #input("Enter the center position as an XYZ file. ")

    x_max = 4 #int(input("Enter the width (x distance) of your image. "))
    y_max = 4 #int(input("Enter the height (y dist) of your image. "))
    edge_buffer = 1 #int(input("Enter the edge buffer distance. "))

    #convert XYZ file (of centers) to list of center objects
    list_of_centers = xyz_to_objects(cpfile)
    
    #make list of all center positions
    positions = []
    for center in list_of_centers:
        position = center.get_location()
        positions.append(position)
    
    # convert data in file into floats and append to a position list
#    with open(cpfile) as f:
#        content = f.readline()
#
#    string = ""
#
#    locations = []
#
#    for i in range(len(content)):
#        if content[i] == " ":
#            locations.append(float(string))
#            string = ""
#        else:
#            string += content[i]
#
#    locations.append(float(string))
#
#    positions = [[""]]
#
#    for i in range(len(locations)):
#        if i % 3 == 0:
#            positions[int(i / 3)] = [locations[i]]
#            positions.append("")
#        else:
#            positions[int(i / 3)].append(locations[i])
#
#    del positions[len(positions) - 1]


    # sort positions for the double finder function

    positions = sorted(positions)

    # Create a Graph of the Input Data
    xypts = []

    for i in range(len(positions)):
        xypts.append([positions[i][0], positions[i][1]])

    # print(xypts)

    points = numpy.array(xypts)
    tri = Delaunay(points)
    #print(len(tri.simplices))

    # print(tri.simplices)

    o_locations = []

    for i in range(len(tri.simplices)):
        midptx1 = 0.50 * (points[tri.simplices][i][0][0] +
                          points[tri.simplices][i][1][0])
        midpty1 = 0.50 * (points[tri.simplices][i][0][1] +
                          points[tri.simplices][i][1][1])
        o_locations.append([midptx1, midpty1, 0])

        midptx2 = (points[tri.simplices][i][1][0] +
                   points[tri.simplices][i][2][0]) / 2.00
        midpty2 = (points[tri.simplices][i][1][1] +
                   points[tri.simplices][i][2][1]) / 2.00
        o_locations.append([midptx2, midpty2, 0])

        midptx3 = (points[tri.simplices][i][2][0] +
                   points[tri.simplices][i][0][0]) / 2.00
        midpty3 = (points[tri.simplices][i][2][1] +
                   points[tri.simplices][i][0][1]) / 2.00
        o_locations.append([midptx3, midpty3, 0])

    print(len(o_locations))
    o_locations.sort
    o_locations = sorted(o_locations)
    # print(o_locations)

    remove = []

    for i in range(len(o_locations) - 1):
        if o_locations[i] == o_locations[i + 1]:
            remove.append(i + 1)

    remove.sort(reverse=True)
    print("Number of O locations: ", len(o_locations))
    # print(remove)

    for i in range(len(remove)):
        del (o_locations[remove[i]])

    print("Number of O locations: ", len(o_locations))


    xOpos = []
    yOpos = []

    for i in range(len(o_locations)):
        xOpos.append(o_locations[i][0])
        yOpos.append(o_locations[i][1])

    # write O positions to an out file
    out = open("OfC Positions 120106_008 Python Output.txt", "w")
    out.write(str(o_locations))
    out.write("nn")

    positions = o_locations

    # find triplets
    triples = o_locator(positions)
#    print("Triples ", triples)

    # find Si positions
    si_locations = []
    for j in range(len(triples)):
        si_locations.append(si_finder(triples[j]))

    # rings = find four(positions, .35)

# --------------------Addition of ring finding------------------------------- #
# assigns nearest 3 adjacent ring to each Si

#    center_objects = []
#    for loc in positions:
#        center = Si_Ring_Classes.ring_center(loc[0], loc[1], loc[2])
#        center_objects.append(center)

    si_objects = []
    for loc in si_locations:
        si = Si_Ring_Classes.Si(loc[0], loc[1], loc[2])
        si.find_rings(list_of_centers, x_max, y_max, edge_buffer)
        print()
        si_objects.append(si)
    
    types = get_stats(si_objects)
    print(types)
    
#    for si in si_objects:
#        print(si.get_location(), end=" ")
#        for ring in si.get_rings():
#            print(ring.get_location(), end=" ")
#        print()

# --------------------------------------------------------------------------- #

    delete = []

    for i in range(len(delete)):
        del si_locations[delete[i]]

    # Plot
    xSipos = []
    ySipos = []

    for i in range(len(si_locations)):
        xSipos.append(si_locations[i][0])
        ySipos.append(si_locations[i][1])

    xOpos = []
    yOpos = []

    for i in range(len(o_locations)):
        xOpos.append(o_locations[i][0])
        yOpos.append(o_locations[i][1])

    plt.triplot(points[:, 0], points[:, 1], tri.simplices.copy())
    plt.plot(points[:, 0], points[:, 1], 'o', color='#2E9AFE')
    plt.scatter(xOpos, yOpos, label='Center Positions', color='#2E9AFE')
    plt.scatter(xOpos, yOpos, label='Oxygen Positions', color='r')
    plt.scatter(xSipos, ySipos, label='Silicon Positions', color='g')

    # write Si positions to an outfile
    out = open("Si Positions Output 170404.txt", "w")
    out.write(str(si_locations))
    out.write("\n")
    plt.xlabel('x (nm)')
    plt.ylabel('y (nm)')
    plt.title('Center Positions')
    plt.legend()
    plt.show()

    # write O positions to an out file
    out = open("OfC Positions 120106_008 Python Output.txt", "w")
    out.write(str(o_locations))
    out.write("nn")


if __name__ == "__main__":
    main()
