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

#import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
    from tkinter import messagebox
else:
    import tkinter as Tk
    from tkinter import messagebox
    
from scipy.spatial import Delaunay
import Si_Ring_Classes
import matplotlib.pyplot as plt


def getNearestNeighbors(base_coords, search_coords, num_neighbors):
    nearest = NearestNeighbors(n_neighbors=num_neighbors, algorithm='ball_tree').fit(base_coords)
    dist, ind = nearest.kneighbors(search_coords)
    return dist, ind

si_locations = numpy.load('Silicon Coords.npy')

si_dist, si_ind = getNearestNeighbors(si_locations, si_locations, 4)


cos_of_angles = []
angles_in_rad = []
final_angle_pairs = []
centers = []
si_dist, si_ind = getNearestNeighbors(si_locations, si_locations, 4)
for i in range(len(si_locations)):
    # si_originvect_dist = si_dist[i][1]
    # si_originvect_ind = [si_locations[si_ind[i][1]][0] - si_locations[si_ind[i][0]][0], 
    #                      si_locations[si_ind[i][1]][1] - si_locations[si_ind[i][0]][1]]
    si_originvect_dist = 1
    si_originvect_ind = [1, 0]
    angle_pair = []
    centersforoneatom = []
    for j in (1, 2, 3):
        # centerangle = []
        si_newvect_dist = si_dist[i][j]
        midpoint_ind = [si_locations[si_ind[i][j]][0] - si_locations[si_ind[i][0]][0] / 2, 
                        si_locations[si_ind[i][j]][1] - si_locations[si_ind[i][0]][1] / 2]
        si_newvect_ind = [si_locations[si_ind[i][j]][0] - si_locations[si_ind[i][0]][0], 
                          si_locations[si_ind[i][j]][1] - si_locations[si_ind[i][0]][1]]
        dot_prod = (si_newvect_ind[0] * si_originvect_ind[0] +
                    si_newvect_ind[1] * si_originvect_ind[1])
        dist_prod = si_originvect_dist * si_newvect_dist
        if si_newvect_ind[1] < 0:
            dot_prod *= -1
        # centerangle.append(centerangle)
        # centerangle.append(dot_prod / dist_prod)
        angle_pair.append(dot_prod / dist_prod)
        centersforoneatom.append(midpoint_ind)
    cos_of_angles.append(angle_pair)
    centers.append(centersforoneatom)
# onlyangles = []
# for i in cos_of_angles:
#     onlyangles.append(i)

for pair in cos_of_angles:
    pairs_in_rad = numpy.arccos(pair)
    for ang in range(len(pair)):
        if pair[ang] < 0:
            pairs_in_rad[ang] *= -1
    angles_in_rad.append(pairs_in_rad)
for pair in angles_in_rad:
    newpair = []
    for angle in pair:
        if angle > (numpy.pi / 2):
            newpair.append((angle * 180 / numpy.pi))
        else:
            newpair.append(angle * 180 / numpy.pi)
    final_angle_pairs.append(newpair)
#print(final_angle_pairs)


x_coords = []
angles = []

for i in range(len(centers)):
    three_centers = centers[i]
    three_angles = final_angle_pairs[i]
    for k in range(len(three_centers)):
        x_coords.append(three_centers[k][0])
        if three_angles[k] < 0:
            three_angles[k] = -180 - three_angles[k]
        angles.append(three_angles[k])
        

plt.scatter(x_coords, angles, s=3)







