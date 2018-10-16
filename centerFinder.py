#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 14:48:58 2018

@author: danwall & thomasmarsh

@editor: rdanilek
"""


import math
import numpy
import matplotlib

matplotlib.use('TkAgg')

from skimage import io
from skimage import feature
from skimage import draw
from skimage import util
from skimage import color
from skimage import morphology
from skimage import filters
from skimage import measure
from skimage import transform
from skimage import exposure
from sklearn.neighbors import NearestNeighbors
from scipy import ndimage as ndi
# import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from scipy.spatial import Delaunay
import Si_Ring_Classes
import matplotlib.pyplot as plt
import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
    from tkinter import messagebox
else:
    import tkinter as Tk
    from tkinter import messagebox



light_centers = False

# 170908_SiO2@RuVII_STM110_det (27.4x11.5).jpg

save_figures = False

# Colors for each of the ring numbers
# DSW edit 5/31/18: made 9 member brown to match PRL
COLORS = [[178, 112, 248], [75, 176, 246], [67, 196, 127],
          [249, 222, 62], [249, 76, 62], [147, 73, 0]]
COLORSHEX = ['#B270F8', '#4BB0F6', '#43C47F', '#F9DE3E', '#F94C3E', '#934900']

scaling_factor = 1

# Sets up the Tk window
root = Tk.Tk()
root.wm_title("Auto Ring Center Finder")

def getFilename():
    """ Runs centerFinder after asking user for inputs """

    # Variables for storing text entered
    entryFilename = Tk.StringVar()
    entryHeight = Tk.StringVar()
    entryWidth = Tk.StringVar()
    entryHoles = Tk.StringVar()
    xyzImport = Tk.IntVar()
    entryXyz = Tk.StringVar()
    prompt = Tk.Label(master=root, text='Enter the filename of the image')
    prompt.pack(side=Tk.TOP)
    filenameEntry = Tk.Entry(master=root, textvariable=entryFilename)
    filenameEntry.pack(side=Tk.TOP)

    def enterFilename():
        # Get Filename from Text Entry
        filename = entryFilename.get()
        xyz_filename = entryXyz.get()
        import_xyz = xyzImport.get()

        try:
            io.imread(filename)
            dimensions = [float(entryWidth.get()), float(entryHeight.get())]
            num_holes = int(entryHoles.get())

            # Destroy the entry UI
            prompt.destroy()
            filenameEntry.destroy()
            continueBtn.destroy()
            widthPrompt.destroy()
            heightPrompt.destroy()
            widthEntry.destroy()
            heightEntry.destroy()
            holesPrompt.destroy()
            holesEntry.destroy()
            xyzPrompt.destroy()
            xyzEntry.destroy()
            importXyzBtn.destroy()

            # Call the rest of the program using given parameters
            centerFinder(filename, dimensions, num_holes, import_xyz,
                         xyz_filename)
        except ValueError:
            messagebox.showerror("Error",
                                 "Invalid dimensions and/or number of holes")
        except FileNotFoundError:
            messagebox.showerror("Error",
                                 "A file could not be found with that filename")

    def importXyz():
        if xyzImport.get():
            xyzPrompt.pack(side=Tk.TOP)
            xyzEntry.pack(side=Tk.TOP)
        else:
            xyzPrompt.pack_forget()
            xyzEntry.pack_forget()

    xyzPrompt = Tk.Label(master=root, text='Enter XYZ Filename')
    xyzEntry = Tk.Entry(master=root, textvariable=entryXyz)

    widthPrompt = Tk.Label(master=root, text='Image Width (nm)')
    widthPrompt.pack(side=Tk.TOP)

    widthEntry = Tk.Entry(master=root, textvariable=entryWidth)
    widthEntry.pack(side=Tk.TOP)

    heightPrompt = Tk.Label(master=root, text='Image Height (nm)')
    heightPrompt.pack(side=Tk.TOP)

    heightEntry = Tk.Entry(master=root, textvariable=entryHeight)
    heightEntry.pack(side=Tk.TOP)

    holesPrompt = Tk.Label(master=root, text='Number of Holes')
    holesPrompt.pack(side=Tk.TOP)

    holesEntry = Tk.Entry(master=root, textvariable=entryHoles)
    holesEntry.pack(side=Tk.TOP)

    importXyzBtn = Tk.Checkbutton(master=root, text='Import XYZ ',
                                  variable=xyzImport, command=importXyz)
    importXyzBtn.pack(side=Tk.TOP)

    continueBtn = Tk.Button(master=root, text='Continue', command=enterFilename)
    continueBtn.pack(side=Tk.TOP)

    # Keeps the window open and listens for events
    while True:
        try:
            root.mainloop()
            break
        # Fixes issue where scrolling accidentally would crash Tkinter
        except UnicodeDecodeError:
            pass

"""AUTOMATIC CENTER FINDING CODE"""


def centerFinder(filename, dimensions, num_holes, import_xyz, xyz_filename):

    # Convert coordinates based on the dimensions of the image
    def pixelsToNm(pixel_coord, nm_dim, im_width, im_height):
        scale = nm_dim[0] / im_width
        return [pixel_coord[0] * scale, pixel_coord[1] * scale]

    def nmToPixels(nm_coord, nm_dim, im_width, im_height):
        scale = im_width / nm_dim[0]
        return [int(nm_coord[0] * scale), int(nm_coord[1] * scale)]

    def pixelDistToNm(dist, nm_dim, im_width, im_height):
        scale = nm_dim[0] / im_width
        return dist * scale

    def importAndScale(filename):
        """Imports and scales the image by a given scaling factor"""
        original = io.imread(filename)
        return transform.rescale(original, scaling_factor, mode='constant',
                                 preserve_range=True).astype('uint8')

    image = importAndScale(filename)

    # Convert image to greyscale
    grey = color.rgb2grey(image)

    # Adjust the brightness of the greyscale to make image more uniform
    grey = exposure.equalize_adapthist(exposure.adjust_gamma(grey))

    # Invert image if centers are dark, else keep same
    if light_centers:
        grey_inv = grey
    else:
        grey_inv = util.invert(grey)

    image_width = len(grey[0])
    image_height = len(grey)

    # Creates instance of STM class
    stm_image = Si_Ring_Classes.STM(filename, (len(grey), len(grey[0])), dimensions, num_holes)

    # Finds the distances and indices of the nearest given number of the base coords
    # to each of the provided search coords
    def getNearestNeighbors(base_coords, search_coords, num_neighbors):

        #num_neighbors = num_neighbors.reshape(1,-1)
        nearest = NearestNeighbors(n_neighbors=num_neighbors,
                                   algorithm='ball_tree').fit(base_coords)

        dist, ind = nearest.kneighbors(search_coords)
        return dist, ind

    def distFromPointToGroup(point, group):
        shortest_dist = math.sqrt(((group[0][0] - point[0])**2) +
                         ((group[0][1] - point[1])**2))
        for coord in group[1:]:
            dist = math.sqrt(((coord[0] - point[0])**2) +
                         ((coord[1] - point[1])**2))
            if dist <= shortest_dist:
                shortest_dist = dist
        return shortest_dist

    stm_image.get_hole_coords(grey)
    holes = stm_image.get_hole_image()

    if not import_xyz:
        #Black out the holes in the greyscale image so no centers will be placed there
        grey_inv = grey_inv * (1 - holes)
        grey_inv_copy = grey_inv

        peak_coord = feature.peak_local_max(grey_inv, min_distance=24,threshold_rel=0.2)

        #Returns a list of the centers of the blobs
        def getBlobCenters(blobs):
            return [[int(blob[1]), int(blob[0])] for blob in blobs]

        #Find blobs in image in order to locate ring centers

        # min_sigma -> increase to find more small rings
        # max_sigma -> increase to find more large rings
        # sigma_ratio -> used in calculation of rings
        # threshold -> decrease to detect less intense rings
        # overlap -> fraction of the blobs that are allowed to overlap with each other

        #blobs = feature.blob_dog(grey_inv, min_sigma=0.03, max_sigma=30,
        #                         sigma_ratio=2.8, threshold=0.8, overlap=0.5)
        blobs = feature.blob_dog(grey_inv, min_sigma=0.07, max_sigma=15,
                                 sigma_ratio=2.8, threshold=0.5, overlap=0.3)
        centers = getBlobCenters(blobs)

        #Find the average distance to closest neighbor
        c_dist, c_ind = getNearestNeighbors(centers, centers, 2)
        avg_closest = numpy.median(c_dist[:][1])
        num_iter = 2

        for i in range(num_iter):
            #Blackout regions around already found ring centers
            #TO DO: Possibly look into scaling this based on average closest
            average_thresh = 1.6 #2.1#1.6 #1.92

            stm_image.plot_circles(grey_inv, centers, avg_closest*average_thresh, 0, False)

            #Find blobs that may be rings that have not been found on first pass
            blobs = numpy.concatenate((blobs, feature.blob_dog(grey_inv,
                                                              min_sigma=0.08,
                                                              max_sigma=20,
                                                              sigma_ratio=2.8,
                                                              threshold=0.8,
                                                              overlap=0.3)))
            centers = getBlobCenters(blobs)

            stm_image.plot_circles(grey_inv_copy, peak_coord,
                                  avg_closest * average_thresh, 0, True)

            peak_coord = numpy.concatenate((peak_coord,
                                            feature.peak_local_max(grey_inv_copy,
                                                                   min_distance=24,
                                                                   threshold_rel=0.2)))

        for k in range(len(peak_coord)):
            temp = peak_coord[k][0]
            peak_coord[k][0] = peak_coord[k][1]
            peak_coord[k][1] = temp

        peak_dist, peak_ind = getNearestNeighbors(centers, peak_coord, 2)

        # Somewhat experimental -> tends to add extra centers but for some images too many
        add_peak_centers = False
        if add_peak_centers:
            peak_thresh = 1.4
            for k in range(len(peak_dist)):
                if peak_dist[k][0] > avg_closest * peak_thresh:
                    centers.append(peak_coord[k])

    else:
        def getCentersFromXyz(xyz_filename):
            centers = []
            with open(xyz_filename, encoding='utf-8') as f:
                file_lines = f.readlines()

            file_lines = [x.strip() for x in file_lines]
            for line in file_lines:
                split_line = line.split(" ")
                nm_coord = [float(split_line[1]) * 10, float(split_line[2]) * 10]
                pixel_coord = nmToPixels(nm_coord, dimensions,
                                         image_width, image_height)
                centers.append(pixel_coord)

            return centers

        centers = getCentersFromXyz(xyz_filename)

    # Find the average distance to closest neighbor
    c_dist, c_ind = getNearestNeighbors(centers, centers, 2)
    average_closest = numpy.median(c_dist[:][1])

    def getNumNeighbors(centers, thresh, average_closest):
        """Returns several lists, including the distance to hole (nm),
        ring sizes, ring center coordinates, and  hole coordinates (nm)"""
        # Get distances and indices of 9 nearest neighbors to every center
        distances, indices = getNearestNeighbors(centers, centers, 10)

        # Gets the distances from every center to the nearest point on the edge of a hole
        if num_holes > 0:
            hole_distances, hole_inds = getNearestNeighbors(stm_image._hole_coords,
                                                            centers, 2)

        hole_dists = [] # List of the distance of each center to nearest hole edge, I think
        ring_size = [] # List of size of each ring
        # center_coord = []  # coordinates of each center (All 3 lists are correlated)

        stm_image._rings = []
        for k in range(len(distances)):
            n_dists = distances[k]

            # Averages the distances to closest 4 centers and multiplies by a
            # threshold to get the max distance for something to be a neighbor
            max_dist = numpy.mean(n_dists[1:5]) * thresh

            #Determines how many of the neighbors are within the max distance
            num_neighbors = 9
            for i in range(4,10):
                if n_dists[i] > max_dist:
                    num_neighbors = i - 1
                    break

            r_full, c_full = draw.circle(centers[k][1], centers[k][0],
                                         max_dist)
            r_bound, c_bound = draw.circle(centers[k][1], centers[k][0],
                                           max_dist, shape=(image_height,
                                                            image_width))

            # Gets the perc of the ring neighbors that are visible in the window
            percent_visible = len(r_full) / len(r_bound)

            # Scales the number of neighbors based on what it should be if all the
            # ring neighbors were visible
            scaled_num_neighbors = int(num_neighbors * percent_visible)

            nearest_centers = []
            for i in range(1, num_neighbors+1):
                nearest_centers.append(centers[indices[k][i]][:])

            #Calculates the centroid of neighboring centers
            x = [p[0] for p in nearest_centers]
            y = [p[1] for p in nearest_centers]
            centroid = (sum(x) / len(nearest_centers), sum(y) / len(nearest_centers))

            exclude_thresh = 1.9

            if num_holes > 0:
                cur_hole_dist = hole_distances[k][1]
            else:
                cur_hole_dist = 2*exclude_thresh*average_closest

            #Gets the coordinates of the circles around the centers
            if (4 <= scaled_num_neighbors <= 9 and percent_visible < 1.2 and
                cur_hole_dist > exclude_thresh*average_closest):
                hole_dists.append(cur_hole_dist)
                ring_size.append(scaled_num_neighbors)
                stm_image._rings.append([(centers[k][0]+centroid[0])/2,(centers[k][1]+centroid[1])/2])

        hole_nm_coords = [] #converts hole coords to nm
        for coord in stm_image._hole_coords:
            hole_nm_coords.append(pixelsToNm(coord, dimensions, image_width, image_height))
            #stm_image._hole_coords.append(pixelsToNm(coord, dimensions, image_width, image_height))

        hole_nm_dists = [] #converts hole_dists to nm
        for dist in hole_dists:
            hole_nm_dists.append(pixelDistToNm(dist, dimensions, image_width, image_height))
            #stm_image._hole_dists.append(pixelDistToNm(dist, dimensions, image_width, image_height))

        return hole_nm_dists, ring_size, hole_nm_coords

    #Threshold for the maximum distance that two centers can be apart to be concidered neighbors
    #thresh = 1.35
    thresh = 1.48
    #thresh = 1.52

    (hole_dist,
     ring_size,
     stm_image._hole_coords) = getNumNeighbors(centers, thresh, average_closest)

    def plotRingCenters(image, ring_size, centers, average_closest):
        for i in range(len(ring_size)):
            # Get circle coordinates for outlines and actual circles for centers
            r_out, c_out = draw.circle(centers[i][1], centers[i][0],
                                       int(average_closest / 3) + 3,
                                       shape=stm_image._im_dim)
            rr, cc = draw.circle(centers[i][1], centers[i][0],
                                 int(average_closest / 3),
                                 shape=stm_image._im_dim)

            # Plot outlines
            image[r_out, c_out] = [0, 0, 0]

            # Assign appropriate COLORS to center coordinates
            image[rr, cc] = COLORS[ring_size[i]-4]

    plotRingCenters(image, ring_size, stm_image._rings, average_closest)

    # Outputs the center data in xyz file format
    def createXyzFile(ring_size):
        text_file = open(filename[:-4]+'xyz.txt', "w")

        for i in range(len(ring_size)):
            nmCoord = pixelsToNm(stm_image._rings[i], dimensions, image_width, image_height)
            text_file.write(str(ring_size[i]),
                            str(nmCoord[0] / 10),
                            str(nmCoord[1] / 10),'0\n')
        text_file.close()

    def distance(position1, position2):
        """finds the distance between two atoms"""
        return math.sqrt((position1[0] - position2[0]) ** 2 +
                         (position1[1] - position2[1]) ** 2 +
                         (position1[2] - position2[2]) ** 2)

    def dists(positions, dist):
        """finds if a triplet could have an Si atom between them"""

        # If there is a triplet close enough to have a Si, return the triplet
        if len(positions) == 3:
            if distance(positions[1], positions[2]) <= dist:
                return positions

        # If there are more than 2 close enough to have a Si between them, find
        # the one that could not be used given the other two
        if len(positions) > 3:
            numbers = []
            for i in range(len(positions)):
                numbers.append(0)
            for i in range(1, len(positions) - 1):
                for j in range(1, len(positions) - i):
                    # If two positions are not close enough, add a counter to both
                    if distance(positions[i], positions[i + j]) > dist:
                        numbers[i] += 1
                        numbers[i + j] += 1
                    # If they are close enough, remove a counter from both
                    else:
                        numbers[i] -= 1
                        numbers[i + j] -= 1

            # Remove the one with the most counters
            del positions[numbers.index(max(numbers))]

            # If close enough, return triplet
            if distance(positions[1], positions[2]) <= dist:
                return positions

        # If they were not enough close to make a triplet, return none
        return[""]

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

    def si_finder(o_pos):
        """finds the position of a Si given a triplet of oxygen"""

        # Characteristic distance
        dist = 1.6 * 0.1

        # Sets up the translation to happen around a basepoint(the first point
        # in the positions)
        trans = [[0, 0, 0], [o_pos[1][0] - o_pos[0][0],
                             o_pos[1][1] - o_pos[0][1],
                             o_pos[1][2] - o_pos[0][2]],
                 [o_pos[2][0] - o_pos[0][0],
                  o_pos[2][1] - o_pos[0][1],
                  o_pos[2][2] - o_pos[0][2]]]

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
        si_pos = [si_pos.item(0) + o_pos[0][0],
                  si_pos.item(1) + o_pos[0][1],
                  si_pos.item(2) + o_pos[0][2]]

        return si_pos

    def o_locator(o_pos):
        """locates all possible triplets"""

        # assumed oxygens are ordered by increasing x values
        # used to collect all the found oxygens close enough to have a single
        # Si between them
        found = [[""]]
        # for each oxygen
        for i in range(len(o_pos)):
            found[i] = [o_pos[i]]
            # for each oxygen with an x position higher than the current
            for j in range(1, len(o_pos) - i):
                # if the x position is less than the possible distance between two
                #  oxygenatoms(variableinclusionradius)
                if abs(o_pos[i][0] - o_pos[i + j][0]) <= 3.45 * 0.1:
                    # if the distance between the two oxygens is less than the
                    #  characteristic distance(variable inclusion radius)
                    if distance(o_pos[i], o_pos[i + j]) <= 3.45 * 0.1:
                        found[i].append(o_pos[i + j])
            found.append([""])

        # removes last appended empty list
        del found[len(found) - 1]

        # remove all those too far apart using dist function (variable
        # inclusion radius)
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

        o_pos = []

        for i in range(len(positions)):
            # center at origin
            pos2 = [positions[i][1][0] - positions[i][0][0],
                    positions[i][1][1] - positions[i][0][1],
                    positions[i][1][2] - positions[i][0][2]]

            # rotate until both points are in the xy plane
            theta = numpy.arctan2(pos2[1], pos2[0])
            phi = numpy.arctan2(pos2[2], pos2[0])
            newx = math.sqrt(math.pow(pos2[0], 2) + math.pow(pos2[2], 2))
            newy = newx * math.tan(theta)

            # find in si position (midpoint between origin and pos 2 in the x-y
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
            o_pos.append([opos[0] + positions[i][0][0],
                          opos[1] + positions[i][0][1],
                          opos[2] + positions[i][0][2]])

        return o_pos

    def centers_to_objects(ring_size, center_list, unit):
        """Converts list of centers to center objects"""
        center_obj_list = []
        for i in range(len(center_list)):
            center = Si_Ring_Classes.ring_center(ring_size[i],
                                                 center_list[i][0],
                                                 center_list[i][1], 0, unit)
            center_obj_list.append(center)
        return center_obj_list

    # stat functions
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
        if len(rings) == 3:
            t1 = int(rings[0].get_type())
            t2 = int(rings[1].get_type())
            t3 = int(rings[2].get_type())
            ordered = order([t1, t2, t3])

            return (str(ordered[0]) + str(ordered[1]) + str(ordered[2]))
        else:
            return "000"

    def get_stats(si_list):
        """ Determines the number of each type of ring triplet [(5, 5, 6),
            (5, 6, 7), etc] and returns a list of tuples and ints containing
            the triplet type and the number found: [(5, 6, 7), 10, (5, 5, 5),
            15, etc]. """
        types = []
        for atom in si_list:
            typ = find_type(atom)

            if typ != '000':
                if typ in types:
                    types[types.index(typ) + 1] += 1
                else:
                    types.append(typ)
                    types.append(1)
        return types

    def getSiOPlot(ring_size):
        x_max = dimensions[0]
        y_max = dimensions[1]
        edge_buffer = 0.1  # Should look into this a bit more

        centers_nm = []
        for coord in stm_image._rings:
            centers_nm.append(pixelsToNm(coord, dimensions,
                                         image_width, image_height))

        # convert XYZ file (of centers) to list of center objects
        list_of_centers = centers_to_objects(ring_size, centers_nm,
                                             'nm')

        # make list of all center positions
        positions = []
        for center in list_of_centers:
            position = center.get_nm_location()
            positions.append(position)

        # sort positions for the double finder function

        positions = sorted(positions)

        # Create a Graph of the Input Data
        xypts = []

        for i in range(len(positions)):
            xypts.append([positions[i][0], positions[i][1]])

        points = numpy.array(xypts)
        tri = Delaunay(points)

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

        o_locations.sort
        o_locations = sorted(o_locations)

        remove = []

        for i in range(len(o_locations) - 1):
            if o_locations[i] == o_locations[i + 1]:
                remove.append(i + 1)

        remove.sort(reverse=True)

        for i in range(len(remove)):
            del (o_locations[remove[i]])

        xOpos = []
        yOpos = []

        for i in range(len(o_locations)):
            xOpos.append(o_locations[i][0])
            yOpos.append(o_locations[i][1])

        positions = o_locations

        # find triplets
        triples = o_locator(positions)

        # find Si positions
        si_locations = []
        for j in range(len(triples)):
            si_locations.append(si_finder(triples[j]))

        # make a list of the Si atoms as objects
        si_objects = []
        for loc in si_locations:
            # construct Si object
            si = Si_Ring_Classes.Si(loc[0], loc[1], loc[2], 'nm')
            # find rings for each object
            si.find_rings(list_of_centers, x_max, y_max, edge_buffer)
            si_objects.append(si)  # add object to list

        types = get_stats(si_objects)  # run stats data and return stat list

        # parse stat list into types and counts
        triplet_types = []
        counts = []
        for i in range(len(types)):
            if i % 2 == 0:
                triplet_types.append(types[i])
            else:
                counts.append(types[i])

        # sort based on most popular type
        triplet_types = [x for _,x in sorted(zip(counts,triplet_types),
                                             reverse=True)]
        counts = sorted(counts, reverse=True)

        # count total
        total_counts = 0
        for count in counts:
            total_counts += count

        plot_image = io.imread(filename)

        xSipos = []
        ySipos = []

        siCoords = []

        for i in range(len(si_locations)):
            siCoords.append(nmToPixels(si_locations[i], dimensions,
                                       image_width, image_height))
            xSipos.append(si_locations[i][0])
            ySipos.append(si_locations[i][1])

        stm_image.plot_circles(plot_image, siCoords, 5, [0, 250, 0], False)

        xOpos = []
        yOpos = []

        oCoords = []

        for i in range(len(o_locations)):
            oCoords.append(nmToPixels(o_locations[i], dimensions,
                                      image_width, image_height))

            xOpos.append(o_locations[i][0])
            yOpos.append(o_locations[i][1])

        # Add oxygens to ring's lists of atoms
        for center in list_of_centers:
            center_position = center.get_nm_location()

            distances, o_inds = getNearestNeighbors(o_locations,
                                                    [center_position],
                                                    center.get_type())

            for o in o_inds[0]:
                center.set_atom(o_locations[o])

        stm_image.plot_circles(plot_image, oCoords, 5, [250, 0, 0], False)

        stm_image.plot_circles(plot_image, stm_image._rings, 5, [0, 0, 250], False)

        return plot_image, o_locations, si_locations, triplet_types, counts, si_objects

    def createWindow(image):
        print()
        print("Launched Auto Ring Finder")
        fig = Figure(figsize=(10, 6), dpi=100)
        ax = fig.add_subplot(111)
        ax.imshow(image)

        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        cid = []

        autoRingToggle = Tk.IntVar()
        autoRingToggle.set(1)

        binSizeTxt = Tk.StringVar()

        """ METHODS FOR MANUALLY EDITING CENTERS"""

        def addcenter(event):
            """Add a center where the user clicked"""
            centers.append([event.xdata, event.ydata])
            replotImage()

        def removecenter(event):
            """Remove the nearest center"""
            min_dist = math.hypot(event.xdata - centers[0][0],
                                  event.ydata - centers[0][1])
            match_ind = 0
            for i in range(len(centers)):
                cur_dist = math.hypot(event.xdata - centers[i][0],
                                      event.ydata - centers[i][1])
                if cur_dist < min_dist:
                    min_dist = cur_dist
                    match_ind = i
            del centers[match_ind]

            replotImage()

        def movecenter(event):
            """Remove the nearest center"""
            min_dist = math.hypot(event.xdata - centers[0][0],
                                  event.ydata - centers[0][1])
            match_ind = 0
            for i in range(len(centers)):
                cur_dist = math.hypot(event.xdata - centers[i][0],
                                      event.ydata - centers[i][1])
                if cur_dist < min_dist:
                    min_dist = cur_dist
                    match_ind = i
            del centers[match_ind]

            # Replaces the center with one in the location clicked
            centers.append([event.xdata, event.ydata]);

            replotImage()

        def replotImage():
            """Replot image"""
            image = io.imread(filename)

            #if autoRingToggle.get():
            (hole_dist,
             ring_size,
             stm_image._hole_coords) = getNumNeighbors(centers, thresh, average_closest)

            plotRingCenters(image, ring_size, stm_image._rings, average_closest)

            ax.imshow(image)
            canvas.draw()

            for i in range(len(cid)):
                canvas.mpl_disconnect(cid[i])

        def ringAdd():
            cid.append(canvas.mpl_connect('button_press_event', addcenter))

        def ringRemove():
            cid.append(canvas.mpl_connect('button_press_event', removecenter))

        def ringMove():
            cid.append(canvas.mpl_connect('button_press_event', movecenter))

        """EDITING BUTTON GUI CODE"""

        addRingBtn = Tk.Button(master=root, text='Add Ring', command=ringAdd)
        addRingBtn.pack(side=Tk.LEFT)

        removeRingBtn = Tk.Button(master=root, text='Remove Ring',
                                  command=ringRemove)
        removeRingBtn.pack(side=Tk.LEFT)

        moveRingBtn = Tk.Button(master=root, text='Move Ring',
                                command=ringMove)
        moveRingBtn.pack(side=Tk.LEFT)

        autoRingChkBtn = Tk.Checkbutton(master=root, text='Autocorrect Rings',
                                        variable=autoRingToggle)
        autoRingChkBtn.pack(side=Tk.LEFT)

        def saveImage():
            print("Saving image...")
            io.imsave(filename[:-4] + 'plotted.jpg', image)
            print("Image Saved")

        saveBtn = Tk.Button(master=root, text='Save Image', command=saveImage)
        saveBtn.pack(side=Tk.RIGHT)

        def xyzFile():
            (hole_dist,
             ring_size,
             stm_image._hole_coords) = getNumNeighbors(centers, thresh,
                                                       average_closest)
            createXyzFile(stm_image._rings, ring_size)

        def doneEditing():
            # Destroy editing buttons
            addRingBtn.destroy()
            removeRingBtn.destroy()
            moveRingBtn.destroy()
            autoRingChkBtn.destroy()
            doneBtn.destroy()

            _, ring_size, _ = getNumNeighbors(centers, thresh,
                                                            average_closest)


            #Plots the silicon and oxygen positions onto image and returns all
            # object needed to construct other plots
            (plot_image,
             o_locations,
             si_locations,
             triplet_types,
             counts,
             si_objects) = getSiOPlot(ring_size)


            # outline hole with light blue
            pix_hole_coords = []
            for coord in stm_image._hole_coords:
                pix_hole_coords.append(nmToPixels(coord,
                                                  dimensions,
                                                  image_width, image_height))
            stm_image.plot_circles(plot_image, pix_hole_coords, 5,
                                   [0, 191, 255], False)

            fig.clf()
            ax = fig.add_subplot(111)
            ax.imshow(plot_image)
            canvas.draw()

            """PLOTTING AND EXPORTING METHODS"""

            def splitRingsIntoBins(bin_size, distances, si_objects):
                """Divides rings into bins based on their distance from a
                hole/line of reference (either a hole or the left side of the
                sample. Takes a bin size, list of distances to hole/line (in
                nm), and a corresponsing list of Si objects.  For each bin,
                this function finds each Si atom within its range, and adds its
                3 associated rings types to the bin.  """
                bin_list = []  # 2-d array of each bin of ring types (types are ints)
                bin_mids = []  # midpoint distance of each ring/bin from edge of hole

                #Get maximum distance from ring to hole or axis (in nm)
                max_dist = float(numpy.max(distances))
                bin_start = 0
                while bin_start < max_dist:
                    # this loop creates the bins

                    # save each midpoint in a list
                    bin_mids.append(bin_start + (bin_size / 2))

                    # put bins in the bin_list
                    # for each bin, set up a list of freqs for each ring type
                    bin_list.append([0, 0, 0, 0, 0, 0])
                    bin_start += bin_size

                for i in range(len(si_objects)):
                    si = si_objects[i] # This is the Si that we are looking at
                    si_dist = distances[i] # Here is its distance from the hole or axis
                    si_rings = si.get_rings() # Here is a list of 3 ring objects for each Si

                    # cycle through bins to find which one to put these 3 ring types in
                    bin_start = 0
                    bin_num = 0
                    while bin_start < max_dist:
                        if bin_start <= si_dist < bin_start + bin_size:
                            # If the distance for the si fits in the bin,
                            for ring in si_rings:
                                # increment each ring type for the si's rings
                                ring_type = ring.get_type()
                                # only add a fraction to the frequency,
                                # ex: 1 Si in a 7-mem ring is 1/7 of the ring
                                bin_list[bin_num][ring_type - 4] += (1 /
                                                                     ring_type)

                        bin_start += bin_size
                        bin_num += 1

                return bin_list, bin_mids

            """Exporting to File"""

            def outputO():
                # write O positions to an out file
                out = open(filename[:-4] + 'Oloc.txt', "w")
                out.write(str(o_locations))
                out.write("\n")

            def outputSi():
                """ write Si positions to an out file. """
                out = open(filename[:-4] + 'Siloc.txt', "w")
                out.write(str(si_locations))
                out.write("\n")

            def tripletBar():
                y_pos = numpy.arange(len(triplet_types))
                plt.title('Triplet Type Histogram')
                plt.bar(y_pos, counts, align='center', alpha=0.5)
                plt.xticks(y_pos, triplet_types, rotation=90)
                plt.show()

            """Drawing Plots"""

            def drawPercPlot(bin_list, bin_mids, use_hole):

                size_perc = [[], [], [], [], [], []]

                for curr_bin in bin_list:
                    curr_bin_tot = sum(curr_bin)  # total num of rings in bin
                    # Divides all ring size tots by the tot num of rings in bin
                    for i in range(6):
                        if curr_bin_tot != 0:
                            size_perc[i].append(curr_bin[i] / curr_bin_tot)
                        else:
                            size_perc[i].append(0)

                # Plots the mid point of the bin vs. percs of each ring size
                for i in range(len(size_perc)):
                    plt.plot(bin_mids, size_perc[i],
                             label=str(i + 4) + ' MR', color=COLORSHEX[i])

                # Adds a legend to the plot
                plt.legend(loc=0, ncol=2)
                plt.title('Ring Size Percentage')
                plt.ylabel('Percentage of Rings')
                if use_hole:
                    plt.xlabel('Distance from Nearest Hole (nm)')
                else:
                    plt.xlabel('Distance from Edge (nm)')
                plt.show()

            def getPercStats(use_hole):
                """See if the bin size is valid, and if so, plot the ring size
                percentages"""
                try:
                    binSize = float(binSizeTxt.get())

                    if use_hole:
                        distancesToHole = []
                        for si in si_objects:
                            dist = distFromPointToGroup(si.get_nm_location(),
                                                        stm_image._hole_coords)
                            distancesToHole.append(dist)

                        bin_list, bin_mids = splitRingsIntoBins(binSize,
                                                                distancesToHole,
                                                                si_objects)
                    else:
                        distancesFromAxis = []
                        for si in si_objects:
                            # create a list of the Si distances from the axis
                            x_coord = si.get_nm_location()[0]
                            distancesFromAxis.append(x_coord)
                        bin_list, bin_mids = splitRingsIntoBins(binSize,
                                                                distancesFromAxis,
                                                                si_objects)

                    drawPercPlot(bin_list, bin_mids, use_hole)

                except ValueError:
                    messagebox.showerror("Error",
                                         "Invalid Bin Size (must be a float)")

            def drawCrystalPlot(bin_list, bin_mids, x_label):
                """Plots the Crystallinity, has a variable horizontal axis
                label to differentiate between hole or axis"""

                # create a list of crystallinities for each bin in bin_list
                crys_list = []
                for curr_bin in bin_list:
                    # divide freq of 6-mem rings by sum of all freqs
                    curr_bin_tot = sum(curr_bin)
                    if curr_bin_tot != 0:
                        # divides by tot # of si
                        crystallinity = curr_bin[2] / (curr_bin_tot)
                    else:
                        crystallinity = 0
                    crys_list.append(crystallinity)

                # plot the list
                plt.plot(bin_mids, crys_list)
                plt.title('Crystallinity')
                plt.xlabel(x_label)
                plt.ylabel('Crystallinity')
                plt.show()

            def getCrystalStats(use_hole):
                """Instructs a function to create bins of ring types associated
                with each Si, then instructs another function to plot them."""
                try:
                    binSize = float(binSizeTxt.get())
                    if use_hole:
                        # Plot in relation to the hole
                        distancesToHole = []
                        for si in si_objects:
                            dist = distFromPointToGroup(si.get_nm_location(),
                                                        stm_image._hole_coords)
                            distancesToHole.append(dist)

                        bin_list, bin_mids = splitRingsIntoBins(binSize,
                                                                distancesToHole,
                                                                si_objects)
                        drawCrystalPlot(bin_list, bin_mids,
                                        'Distance from Hole (nm)')
                    else:
                        # Plot in relation to the axis
                        distancesFromAxis = []
                        for si in si_objects:
                            # create a list of every Si distance from the axis
                            x_coord = si.get_nm_location()[0]
                            distancesFromAxis.append(x_coord)

                        bin_list, bin_mids = splitRingsIntoBins(binSize,
                                                                distancesFromAxis,
                                                                si_objects)

                        drawCrystalPlot(bin_list, bin_mids,
                                        'Distance from Edge (nm)')
                except ValueError:
                    messagebox.showerror("Error",
                                       "Invalid Bin Size (must be a float)")

            def ddoPlot(si_locations, hole_coords, use_hole_dist):
                angle_coords = []
                angles = []
                x_coords = []

                si_dist, si_ind = getNearestNeighbors(si_locations,
                                                      si_locations, 4)

                for i in range(len(si_locations)):
                    start_coord = si_locations[si_ind[i][0]]
                    for j in (1, 2, 3):
                        end_coord = si_locations[si_ind[i][j]]
                        ddo_vector = [end_coord[0] - start_coord[0],
                                      end_coord[1] - start_coord[1]]
                        slope = ddo_vector[1] / ddo_vector[0]
                        midpoint_coord = [(start_coord[0] + end_coord[0]) / 2,
                                          (start_coord[1] + end_coord[1]) / 2]
                        rad_angle = math.atan(slope)
                        degree_angle = (rad_angle * 180 / numpy.pi)
                        angles.append(degree_angle)
                        angle_coords.append(midpoint_coord)
                        x_coords.append(midpoint_coord[0])

                if len(hole_coords) > 0:
                    hole_distances, hole_inds = getNearestNeighbors(hole_coords, angle_coords, 1)
                    angle_dists = []
                    for dist in hole_distances:
                        angle_dists.append(dist[0])

                if use_hole_dist:
                    plt.scatter(angle_dists, angles, s=2, marker='_')
                else:
                    plt.scatter(x_coords, angles, s=2, marker='_')

                plt.title('DDO Plot')
                if use_hole_dist:
                    plt.xlabel('Distance from Hole (nm)')
                else:
                    plt.xlabel('x coordinate (nm)')

                plt.ylabel('Angle (degrees)')
                plt.show()

            def percPlotHole():
                getPercStats(True)

            def percPlotX():
                getPercStats(False)

            def plotCrystalHoles():
                getCrystalStats(True)

            def plotCrystalInX():
                getCrystalStats(False)

            def plotDDOHoles():
                ddoPlot(si_locations, stm_image._hole_coords, True)

            def plotDDOinX():
                ddoPlot(si_locations, stm_image._hole_coords, False)

            """PLOTTING AND EXPORTING DROPDOWN MENU CODE"""

            menubar = Tk.Menu(root)

            plotMenu = Tk.Menu(menubar, tearoff=0)
            plotMenu.add_command(label="Ring Size Percentage from Hole",
                                 command=percPlotHole)
            plotMenu.add_command(label="Ring Size Percentage in X",
                                 command=percPlotX)
            plotMenu.add_command(label="Crystallinity from Hole",
                                 command=plotCrystalHoles)
            plotMenu.add_command(label="Crystallinity in X",
                                 command=plotCrystalInX)
            plotMenu.add_command(label="Triplet Type Bar",
                                 command=tripletBar)
            plotMenu.add_command(label="DDO Plot from Hole",
                                 command=plotDDOHoles)
            plotMenu.add_command(label="DDO Plot in X", command=plotDDOinX)
            menubar.add_cascade(label="Plot", menu=plotMenu)

            # Creates the export dropdown menu
            exportMenu = Tk.Menu(menubar, tearoff=0)
            exportMenu.add_command(label="XYZ File", command=xyzFile)
            exportMenu.add_command(label="Oxygen Locations", command=outputO)
            exportMenu.add_command(label="Silicon Locations", command=outputSi)
            menubar.add_cascade(label="Export", menu=exportMenu)

            root.config(menu=menubar)

            # Label for bin size text entry
            binSizeLabel = Tk.Label(master=root, text='Bin Size (nm)')
            binSizeLabel.pack(side=Tk.LEFT)

            # Bin size text entry field
            binSizeEntry = Tk.Entry(master=root, textvariable=binSizeTxt)
            binSizeEntry.pack(side=Tk.LEFT)

        # Button to finish editing and move to exporting / plotting graphs
        doneBtn = Tk.Button(master=root, text='Done Editing',
                            command=doneEditing)
        doneBtn.pack(side=Tk.RIGHT)

    createWindow(image)


getFilename()
