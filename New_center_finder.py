#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 16:16:52 2018

@author: catherineryczek
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


""" Globals """
light_centers = False
save_figures = False
# Colors for each of the ring numbers
# DSW edit 5/31/18: made 9 member brown to match PRL
COLORS = [[178, 112, 248], [75, 176, 246], [67, 196, 127],
          [249, 222, 62], [249, 76, 62], [147, 73, 0]]
scaling_factor = 1

# Sets up the Tk window
root = Tk.Tk()
root.wm_title("Auto Ring Center Finder")


class UserInterface():
    """ A starting class that gathers the preliminary information inputted in
        the user interface and deals with interactions with the user """
    
    def __init__(self):
        """ Constructor """
        
        self._filename = Tk.StringVar()
        self._width = Tk.StringVar()
        self._height = Tk.StringVar()
        self._num_holes = Tk.StringVar()
        self._xyz_file = Tk.StringVar()
        
        # Image filename
        self._prompt = Tk.Label(master=root, text='Enter the filename of the image')
        self._prompt.pack(side=Tk.TOP)
        self._filename_entry = Tk.Entry(master=root, textvariable=self._filename)
        self._filename_entry.pack(side=Tk.TOP)
    
        # Width in nm
        self._width_prompt = Tk.Label(master=root, text='Image Width (nm)')
        self._width_prompt.pack(side=Tk.TOP)
        self._width_entry = Tk.Entry(master=root, textvariable=self._width)
        self._width_entry.pack(side=Tk.TOP)
    
        # Height in nm
        self._height_prompt = Tk.Label(master=root, text='Image Height (nm)')
        self._height_prompt.pack(side=Tk.TOP)
        self._height_entry = Tk.Entry(master=root, textvariable=self._height)
        self._height_entry.pack(side=Tk.TOP)
    
        # Number of holes
        self._holes_prompt = Tk.Label(master=root, text='Number of Holes')
        self._holes_prompt.pack(side=Tk.TOP)
        self._holes_entry = Tk.Entry(master=root, textvariable=self._num_holes)
        self._holes_entry.pack(side=Tk.TOP)
    
        # XYZ file
        self._xyz_prompt = Tk.Label(master=root, text='Enter XYZ Filename')
        self._xyz_entry = Tk.Entry(master=root, textvariable=self._xyz_file)

        # XYZ checkbox button
        self._xyz_import = Tk.IntVar()
        self._import_xyz_btn = Tk.Checkbutton(master=root, text='Import XYZ ',
                                      variable=self._xyz_import,
                                      command=self.importXyz)
        self._import_xyz_btn.pack(side=Tk.TOP)
    
        # Continue button
        self._continue_btn = Tk.Button(master=root, text='Continue', command=self.enterFilename)
        self._continue_btn.pack(side=Tk.TOP)
    
        # Keeps the window open and listens for events
        while True:
            try:
                root.mainloop()
                break
            # Fixes issue where scrolling accidentally would crash Tkinter
            except UnicodeDecodeError:
                pass
        

    def importXyz(self):
        """ Import XYZ data if available """
                
        if self._xyz_file.get():
            self._xyz_prompt.pack(side=Tk.TOP)
            self._xyz_entry.pack(side=Tk.TOP)
        else:
            self._xyz_prompt.pack_forget()
            self._xyz_entry.pack_forget()


    def enterFilename(self):
        """ Get Filename from text entry """
                
        self._filename = self._filename.get()
    
        try:
            io.imread(self._filename)
    
            # Destroy the entry UI
            self._prompt.destroy()
            self._filename_entry.destroy()
            self._continue_btn.destroy()
            self._width_prompt.destroy()
            self._height_prompt.destroy()
            self._width_entry.destroy()
            self._height_entry.destroy()
            self._holes_prompt.destroy()
            self._holes_entry.destroy()
            self._xyz_prompt.destroy()
            self._xyz_entry.destroy()
            self._import_xyz_btn.destroy()
        
            # Pass information on to the main class
            self._image = Image(self._filename, self._width.get(), self._height.get(),
                                self._num_holes.get(), self._xyz_import.get(),
                                self._xyz_file.get(), self)

        except ValueError:
            messagebox.showerror("Error", "Invalid dimensions and/or number of holes")
        except FileNotFoundError:
            messagebox.showerror("Error", "A file could not be found with that filename")
                
        
    def createWindow(self, image, centers):
        """ Create and manage the secondary window allowing for manual image
            editing and advancement to the final display """
        
        self._image_file = image
        self._centers = centers
        
        print()
        print("Launched Auto Ring Finder")
        fig = Figure(figsize=(10, 6), dpi=100)
        self._ax = fig.add_subplot(111)
        self._ax.imshow(self._image_file)

        self._canvas = FigureCanvasTkAgg(fig, master=root)
        self._canvas.draw()
        self._canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(self._canvas, root)
        toolbar.update()
        self._canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        self._cid = []

        autoRingToggle = Tk.IntVar()
        autoRingToggle.set(1)

        self._binSizeTxt = Tk.StringVar()
        
        # EDITING BUTTON GUI CODE

        self._addRingBtn = Tk.Button(master=root, text='Add Ring', command=self.ringAdd)
        self._addRingBtn.pack(side=Tk.LEFT)

        self._removeRingBtn = Tk.Button(master=root, text='Remove Ring',
                                        command=self.ringRemove)
        self._removeRingBtn.pack(side=Tk.LEFT)

        self._moveRingBtn = Tk.Button(master=root, text='Move Ring',
                                      command=self.ringMove)
        self._moveRingBtn.pack(side=Tk.LEFT)

        self._autoRingChkBtn = Tk.Checkbutton(master=root, text='Autocorrect Rings',
                                              variable=autoRingToggle)
        self._autoRingChkBtn.pack(side=Tk.LEFT)
        
        saveBtn = Tk.Button(master=root, text='Save Image', command=self.saveImage)
        saveBtn.pack(side=Tk.RIGHT)
        
        # Button to finish editing and move to exporting / plotting graphs
        self._doneBtn = Tk.Button(master=root, text='Done Editing',
                                  command=self.doneEditing)
        self._doneBtn.pack(side=Tk.RIGHT)
    
    
    """ METHODS FOR MANUALLY EDITING CENTERS """
    
    def replotImage(self):
        """ Replot the image after editing """
        
        # Get a clean copy of the image
        self._image_file = io.imread(self._filename)
        
        self._image.replotImage(self._image_file, self._centers)

        self._ax.imshow(self._image_file)
        self._canvas.draw()

        for i in range(len(self._cid)):
            self._canvas.mpl_disconnect(self._cid[i])
    
    
    def ringAdd(self):
        self._cid.append(self._canvas.mpl_connect('button_press_event', self.addCenter))
    
    
    def addCenter(self, event):
        """ Add a center where the user clicked """
        
        self._centers.append(Center([event.xdata, event.ydata]))
        self.replotImage()
    
    
    def ringRemove(self):
        self._cid.append(self._canvas.mpl_connect('button_press_event', self.removeCenter))
    
    
    def removeCenter(self, event, move_center=False):
        """ Remove the nearest center to location of click event """
        
        # Identify the nearest center to location clicked
        min_dist = math.hypot(event.xdata - self._centers[0].getCoords()[0],
                              event.ydata - self._centers[0].getCoords()[1])
        match_ind = 0
        for i in range(len(self._centers)):
            cur_dist = math.hypot(event.xdata - self._centers[i].getCoords()[0],
                                  event.ydata - self._centers[i].getCoords()[1])
            if cur_dist < min_dist:
                min_dist = cur_dist
                match_ind = i

        del self._centers[match_ind]

        # Replot the image unless a new center needs to be added
        if not move_center:
            self.replotImage()
    
    
    def ringMove(self):
        self._cid.append(self._canvas.mpl_connect('button_press_event', self.moveCenter))
    
    
    def moveCenter(self, event):
        """ Remove the nearest center and add a new center """
        
        # Remove the nearest center
        self.removeCenter(event, True)
        
        # Add a center in the location clicked
        self.addCenter(event)
    
    
    def saveImage(self):
        """ Saves the image """
        
        print("Saving image...")
        io.imsave(self._filename[:-4] + 'plotted.jpg', self._image_file)
        print("Image Saved")
    
    
    def doneEditing(self):
        """ Run the main program after done editing centers """
        
        # Destroy editing buttons
        self._addRingBtn.destroy()
        self._removeRingBtn.destroy()
        self._moveRingBtn.destroy()
        self._autoRingChkBtn.destroy()
        self._doneBtn.destroy()
        
        # Plot silicon and oxygen positions onto image
        self._image.getSiOPlot()



class Image():
    """ An overarching class to describe the STM image. Includes information
        like filename, image dimensions (pixels), sample dimensions (nm),
        scale, number of holes, and coordinates of those holes."""

    def __init__(self, filename, width, height, num_holes, xyz_import, xyz_file,
                 user_interface):
        """ Constructor """

        self._filename = filename
        self._sample_dim = [float(width), float(height)]  # (nm)
        self._num_holes = int(num_holes)
        self._xyz_import = xyz_import
        self._xyz_file = xyz_file
        self._user_interface = user_interface
        
        self._image = self.importAndScale()  # The actual image
        
        self._grey = ''  # Greyscale version of the image
        self._grey_inv = ''  # Inverted version of grey
        
        self.grey()
        
        self._im_dim = (len(self._grey), len(self._grey[0])) # (image height, image width) (pixels)
        self._scale = self._im_dim[1] / self._sample_dim[0]  # Ratio pixels/nm

        # Initialize the hole class
        self._holes = Holes(self._num_holes, self._grey, self._im_dim, self)

        # Initialize the rings class
        self._rings = Rings(self._num_holes, self._holes, self._im_dim,
                            self._xyz_import, self._xyz_file, self._grey_inv, self)
        
        # Return to the user interface for more interactions
        self._user_interface.createWindow(self._image, self._rings.returnCenters())
        

    def nmToPixels(self, nm_coord):
        """ Converts a set of nm coordinates into pixel coordinates """
        return [int(nm_coord[0] * self._scale), int(nm_coord[1] * self._scale)]
    
    
    def pixelsToNm(self, pixel_coord):
        """ Converts a set of pixel coordinates into nm coordinates """
        return [pixel_coord[0] / self._scale, pixel_coord[1] / self._scale]
    
    
    def pixelDistToNm(self, dist):
        """Converts a distance in pixels into a distance in nm """
        return dist / self._scale
    
    
    def distance(self, position1, position2):
        """ Returns the distance between two atoms """
        return math.sqrt((position1[0] - position2[0]) ** 2 +
                         (position1[1] - position2[1]) ** 2 +
                         (position1[2] - position2[2]) ** 2)
    
    
    def returnHoleImage(self):
        return self._holes.returnHoleImage()
    
    
    def importAndScale(self):
        """ Imports and scales the image by a given scaling factor """
        
        original = io.imread(self._filename)
        return transform.rescale(original, scaling_factor, mode='constant',
                                 preserve_range=True).astype('uint8')


    def grey(self):
        """ Convert the imported image to greyscale and adjust exposure """
        
        # Create a greyscale version of the image
        self._grey = color.rgb2grey(self._image)
        
        # Adjust the brightness of the greyscale to make image more uniform
        self._grey = exposure.equalize_adapthist(exposure.adjust_gamma(self._grey))
    
        # Refering the the set global variable, invert image if centers are
        # dark, else keep same
        if light_centers:
            self._grey_inv = self._grey
        else:
            self._grey_inv = util.invert(self._grey)
    
    
    def getNearestNeighbors(self, base_coords, search_coords, num_neighbors):
        """ Finds the distances and indices of the nearest given number of base
            coords to each of the provided search coords. """

        nearest = NearestNeighbors(n_neighbors=num_neighbors,
                                   algorithm='ball_tree').fit(base_coords)
        dist, ind = nearest.kneighbors(search_coords)

        return dist, ind
    
    
    def plotCirclesOnImage(self, image, lst, radius, color, reverse=False):
        """ Plots circles of a given radius and color on a given image at given coordinates """
        if not reverse:
            for center in lst:
                coord = center.getCoords()
                rr, cc = draw.circle(coord[1], coord[0], radius, shape=self._im_dim)
                image[rr, cc] = color
        else:
            # Same as above but coord reversed - deals with peak_cord not center objects
             for coord in lst:
                rr, cc = draw.circle(coord[0], coord[1], radius, shape=self._im_dim)
                image[rr, cc] = color
    
    
    def plotRingCenters(self, center_objects, average_closest):
        """ Plots circles on image of the correct color for the ring size """
        
        for i in range(len(center_objects)):
            # Get circle coordinates for outlines and actual circles for centers
            r_out, c_out = draw.circle(center_objects[i].getRingCenterCoord()[1],
                                       center_objects[i].getRingCenterCoord()[0],
                                       int(average_closest / 3) + 3,
                                       shape=self._im_dim)
            rr, cc = draw.circle(center_objects[i].getRingCenterCoord()[1],
                                 center_objects[i].getRingCenterCoord()[0],
                                 int(average_closest / 3),
                                 shape=self._im_dim)

            # Plot outlines
            self._image[r_out, c_out] = [0, 0, 0]

            # Assign appropriate colors to center coordinates
            self._image[rr, cc] = COLORS[center_objects[i].getRingSize() - 4]
    
    
    def replotImage(self, image, centers):
        """ Replots the image after editing """
        
        self._image = image
        self._rings.updateCenters(centers)
        
        # Rerun ring determination with new centers and replot the new image
        self._rings.centerInfo()
    
    
    def getSiOPlot(self):
        """ The main program after done editing centers """
        
        self._rings.getSiOPlot()
        
    

class Holes():
    """ A class to keep track of hole related information """
    
    def __init__(self, num_holes, grey, image_dim, image_class):
        """ Constructor """
        
        self._num_holes = num_holes
        self._grey = grey  # Grey scale image
        self._im_dim = image_dim  # (image height, image width) (pixels)
        self._image = image_class
        
        self._hole_coords = self.getHoleCoords()  # In pixels
        self._hole_image = self.getHoleImage()  # Mask image of holes
        
        self._hole_nm_coords = []  # In nm
        for coord in self._hole_coords:
            self._hole_nm_coords.append(self._image.pixelsToNm(coord))
        
    
    def getCoords(self):
        return self._hole_coords
    
    
    def returnHoleImage(self):
        return self._hole_image
    
    
    def getHoleCoords(self):
        """ Get coordinates of borders of holes from greyscale image """

        if self._num_holes > 0:
            # Threshold greyscale image to create a binary image
            im_thresh = filters.threshold_minimum(self._grey)
            binary = self._grey > im_thresh

            # Erode the binary image and then subtract from original to get borders
            erosion = numpy.pad(morphology.binary_erosion(binary)[2:-2,2:-2],2,'maximum')
            borders = binary ^ erosion
            borders = numpy.pad(borders[1:-1,1:-1],1,'edge')

            # Label the regions in the image
            label_image = measure.label(borders)
            regions = measure.regionprops(label_image)

            # Put all of the areas of regions in border image into list
            areas = []
            for region in regions:
                areas.append(region.area)

            # Get the coordinates of the largest num_holes number of border regions
            hole_coords = []
            for i in range(self._num_holes):
                max_ind = numpy.argmax(areas)
                coords = regions[max_ind].coords
                for coord in coords:
                    hole_coords.append([coord[1], coord[0]])

                del regions[max_ind]
                del areas[max_ind]

            return hole_coords
        
        else:
            return []
    
    
    def getHoleImage(self):
        """ Takes the coordinates of the borders of the holes and fills them in
            to get a mask image of the holes """

        hole_image = numpy.zeros(self._im_dim)
        
        for coord in self._hole_coords:
            hole_image[coord[1], coord[0]] = 1

        return ndi.binary_fill_holes(morphology.binary_closing(hole_image))



class Rings():
    """ A class dealing with interpreting the rings within the STM image """

    def __init__(self, num_holes, hole_object, image_dimensions, xyz_import, xyz_file,
                 grey_inv, image_class):
        """ Constructor """
        
        self._num_holes = num_holes
        self._holes = hole_object
        self._im_dim = image_dimensions
        self._xyz_import = xyz_import
        self._xyz_file = xyz_file
        self._grey_inv = grey_inv
        self._image = image_class
        
        self._centers = []  # List of center objects
        self._oxygens = []  # List of all oxygen objects
        self._o_trip_locs = []  # List of oxygen triplet locations
        
        # If not importing an XYZ file, identify Si centers on greyscale image
        if not self._xyz_import:
            self._peak_coord = self.blackOutHoles(self._image.returnHoleImage())
            self.getCentersFromImage()
        # Else, retrieve them from XYZ file
        else:
            self.getCentersFromXYZ()
        
        # Threshold for the maximum distance two centers can be apart to be neighbors
        self._neighbor_thresh = 1.48
        
        # Determine each center's distance from the hole(s), ring size, and
        # ring center coordinate and keep in list attributes. Then plot rings
        self.centerInfo()
        
    
    def returnCenters(self):
        return self._centers
    
    
    def updateCenters(self, new_centers):
        """ Redefine the list of centers """
        self._centers = new_centers
    
    
    def blackOutHoles(self, holes):
        """ Black out the holes in the greyscale image so no centers will be
            placed there """

        self._grey_inv = self._grey_inv * (1 - holes)
        self._grey_inv_copy = self._grey_inv

        return feature.peak_local_max(self._grey_inv, min_distance=24, threshold_rel=0.2)
        
    
    def getCentersFromImage(self):
        """ Returns a list of the centers of the rings from the greyscale image"""
        
        blobs = feature.blob_dog(self._grey_inv, min_sigma=0.07, max_sigma=15,
                                 sigma_ratio=2.8, threshold=0.5, overlap=0.3)

        self.getBlobCenters(blobs)
        
        # Find the average distance to closest neighbor
        centers = []
        for center in self._centers:
            centers.append(center.getCoords())
        c_dist, c_ind = self._image.getNearestNeighbors(centers, centers, 2)
        avg_closest = numpy.median(c_dist[:][1])
        num_iter = 2
        
        for i in range(num_iter):
            #Blackout regions around already found ring centers
            #TO DO: Possibly look into scaling this based on average closest
            average_thresh = 1.6

            self._image.plotCirclesOnImage(self._grey_inv, self._centers, avg_closest*average_thresh, 0)

            #Find blobs that may be rings that have not been found on first pass
            blobs = numpy.concatenate((blobs, feature.blob_dog(self._grey_inv, min_sigma=0.08,
                                                              max_sigma=20, sigma_ratio=2.8,
                                                              threshold=0.8, overlap=0.3)))
            self.getBlobCenters(blobs)

            self._image.plotCirclesOnImage(self._grey_inv_copy, self._peak_coord,
                                           avg_closest * average_thresh, 0, True)

            self._peak_coord = numpy.concatenate((self._peak_coord,
                                                  feature.peak_local_max(self._grey_inv_copy, 
                                                     min_distance=24, threshold_rel=0.2)))
                    
    
    def getBlobCenters(self, blobs):
        """ Returns a list of the centers of the blobs """
        
        # Find blobs in image in order to locate ring centers:
        # min_sigma -> increase to find more small rings
        # max_sigma -> increase to find more large rings
        # sigma_ratio -> used in calculation of rings
        # threshold -> decrease to detect less intense rings
        # overlap -> fraction of the blobs that are allowed to overlap with each other

        self._centers = []  # Start with a clean list
        for blob in blobs:
            x_coord = int(blob[1])
            y_coord = int(blob[0])
            self._centers.append(Center((x_coord, y_coord), self._image))
    
    
    def getCentersFromXyz(self):
        """ Read the XYZ file to find the centers noted within """
        
        with open(self._xyz_file, encoding='utf-8') as f:
            file_lines = f.readlines()

        file_lines = [x.strip() for x in file_lines]
        for line in file_lines:
            split_line = line.split(" ")
            nm_coord = [float(split_line[1]), float(split_line[2])]
            pixel_coord = self._image.nmToPixels(nm_coord)
            self._centers.append(Center(pixel_coord), self._image)
    
    
    def centerInfo(self):
        """ Assigns several attributes to each center, including the distance
            to hole (nm), ring sizes, and ring center coordinates with respect
            to each center """
        
        exclude_thresh = 1.9
        
        # Make a current list of center locations
        centers = []
        for center in self._centers:
            centers.append(center.getCoords())
        
        # Get distances and indices of 9 nearest neighbors to every center
        distances, indices = self._image.getNearestNeighbors(centers, centers, 10)
        
        # Find the average distance to closest neighbor
        c_dist, c_ind = self._image.getNearestNeighbors(centers, centers, 2)
        average_closest = numpy.median(c_dist[:][1])
        
        # Gets the distances from every center to the nearest point on the edge of a hole
        if self._num_holes > 0:
            hole_distances, hole_inds = self._image.getNearestNeighbors(self._holes.getCoords(),
                                                                        centers, 2)

        for k in range(len(distances)):
            num_neighbors, percent_visible = self.findRingSizes(k, distances[k])
            
            # Scales num neighbors to be as if all ring neighbors were visible
            scaled_num_neighbors = int(num_neighbors * percent_visible)
            
            # Coordinate of center's ring center
            ring_center_coord = self.findRingCenterCoord(k, num_neighbors, indices)
            
            # Distance from hole(s) to center
            if self._num_holes > 0:
                hole_dist = hole_distances[k][1]
            else:
                hole_dist = 2 * self._exclude_thresh * average_closest
                
            # Convert the pixel distance into nm
            hole_nm_dist = self._image.pixelDistToNm(hole_dist)

            # If all numbers fits the needed criteria
            if 4 <= scaled_num_neighbors <= 9 and percent_visible < 1.2 and\
            hole_dist > exclude_thresh * average_closest:
                self._centers[k].assignRingSize(scaled_num_neighbors)
                self._centers[k].assignRingCenterCoord(ring_center_coord)
                self._centers[k].assignHoleDist(hole_nm_dist)
        
        # Remove centers that did not meet the critera
        new_centers = []
        for center in self._centers:
            if center.getRingCenterCoord():
                new_centers.append(center)
        self._centers = new_centers

        # Plot the rings on the image
        self._image.plotRingCenters(self._centers, average_closest)
    
    
    def findRingSizes(self, k, n_dists):
        """ Determine the ring size by looking at the number of nearest
            neighbors it has """

        # Averages the distances to closest 4 centers and multiplies by a
        # threshold to get the max distance for something to be a neighbor
        max_dist = numpy.mean(n_dists[1:5]) * self._neighbor_thresh

        # Determines how many of the neighbors are within the max distance
        num_neighbors = 9
        for i in range(4,10):
            if n_dists[i] > max_dist:
                num_neighbors = i - 1
                break

        r_full, c_full = draw.circle(self._centers[k].getCoords()[1],
                                     self._centers[k].getCoords()[0], max_dist)
        r_bound, c_bound = draw.circle(self._centers[k].getCoords()[1],
                                       self._centers[k].getCoords()[0],
                                       max_dist, shape=self._im_dim)

        # Percentage of the ring neighbors that are visible in the window
        percent_visible = len(r_full) / len(r_bound)
            
        return num_neighbors, percent_visible
    
    
    def findRingCenterCoord(self, k, num_neighbors, indices):
        """ Determine the coordinates of each center's ring """
        
        nearest_centers = []
        for i in range(1, num_neighbors+1):
            nearest_centers.append(self._centers[indices[k][i]].getCoords()[:])

        #Calculates the centroid of neighboring centers
        x = [p[0] for p in nearest_centers]
        y = [p[1] for p in nearest_centers]
        centroid = (sum(x) / len(nearest_centers), sum(y) / len(nearest_centers))
        
        return [(self._centers[k].getCoords()[0] + centroid[0]) / 2,
                (self._centers[k].getCoords()[1] + centroid[1]) / 2]
    
    
    def getSiOPlot(self):
        """ Obtain information necessary for locating Si and O atoms """
        
        self.makeOxygens()
        self.oLocator()
        self.siLocator()
    
    
    def makeOxygens(self):
        """ Make a list of oxygen objects """
        
        # Create a sorted list of ring center coordinates in nm
        ring_center_nm_coords = []
        for center in self._centers:
            ring_center_nm_coords.append([center.getRingCenterNmCoord()[0],
                                          center.getRingCenterNmCoord()[1]])

        ring_center_nm_coords = sorted(ring_center_nm_coords)
        
        points = numpy.array(ring_center_nm_coords)
        tri = Delaunay(points)
        
        # Make a list of oxygen locations
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
        
        o_locations = sorted(o_locations)

        # Check and remove any duplicates
        new_o_locations = []
        for i in range(len(o_locations) - 1):
            if o_locations[i] != o_locations[i + 1]:
                new_o_locations.append(o_locations[i])
        
        for location in new_o_locations:
            self._oxygens.apppend(Oxygen(location))
        
        
    def oLocator(self):
        """ Locates all possible triplets assumeing oxygens are ordered by
            increasing x values. Used to mark all the found oxygens close
            enough to have a single Si between them """
        
        inclusion_radius = 0.345  # Max possible distance between two oxygen atoms (nm)
        
        for i in range(len(self._oxygens)):
            
            # For each oxygen with an x position higher than the current
            for j in range(1, len(self._oxygens) - i):
                
                # Preliminary calculation - is deltat x less than inclusion radius
                if abs(self._oxygens[i].getCoords[0] - self._oxygens[i + j].getCoords[0])\
                        <= inclusion_radius:
                    
                    # Is actual dist between the two oxygens less than inclusion radius
                    if self._image.distance(self._oxygens[i].getCoords[0], self._oxygens[i + j].getCoords[0]())\
                            <= inclusion_radius:
                        
                        # Add the acceptably close oxygen to the oxygen's triplet list
                        self._oxygens[i].addTriplet(self._oxygens[i + j])
            
            # Check each oxygen's possible triplets to keep only the best
            self.tripleChecker(self._oxygens[i], inclusion_radius)
        
            if self._oxygens[i].getTripletPositions():
                self._o_trip_locs.append(self._oxygens[i].getTripletPositions())
    
    
    def tripleChecker(self, oxygen, inclusion_radius):
        """ Determine the best oxygen triplet to have an Si atom """
        
        # List of possible triplet oxygen atom positions
        positions = []
        for o_atom in oxygen.getTriplet():
            positions.append(o_atom.getCoords())
        
        # If there is a triplet close enough to have a Si, return the triplet
        if len(positions) == 3:
            if self._image.distance(positions[1], positions[2]) <= inclusion_radius:
                oxygen.setTripletPositions(positions)
                return
            
        # If there are more than 2 close enough to have a Si between them, find
        # the one that could not be used given the other two
        if len(positions) > 3:
            numbers = []
            for i in range(len(positions)):
                numbers.append(0)
            for i in range(1, len(positions) - 1):
                for j in range(1, len(positions) - i):
                    # If two positions are not close enough, add a counter to both
                    if self._image.distance(positions[i], positions[i + j]) > inclusion_radius:
                        numbers[i] += 1
                        numbers[i + j] += 1
                    # If they are close enough, remove a counter from both
                    else:
                        numbers[i] -= 1
                        numbers[i + j] -= 1
    
            # Remove the one with the most counters
            del positions[numbers.index(max(numbers))]
                
            # If close enough, return triplet
            if self._image.distance(positions[1], positions[2]) <= inclusion_radius:
                oxygen.setTripletPositions(positions)
                return
    
        # If they were not enough close to make a triplet, return none
        oxygen.setTripletPositions('')
    
    
    def siLocator(self):
        """ Create a list of Si locations and create Si objects in Si_Ring_Classes.py """
        
        si_locations = []
        for triplet in self._o_trip_locs:
            si_locations.append(self.siFinder(triplet))
        
        # TO DO: Look into this more
        # edge_buffer = 0.1
    
    
    def siFinder(self, triplet):
        """ Find the position of a Si given a triplet of oxygen """

        dist = 0.16   # Characteristic distance (nm)

        # Sets up the translation to happen around a basepoint (the first point
        # in the positions)
        trans = [[0, 0, 0], [triplet[1][0] - triplet[0][0],
                             triplet[1][1] - triplet[0][1],
                             triplet[1][2] - triplet[0][2]],
                 [triplet[2][0] - triplet[0][0],
                  triplet[2][1] - triplet[0][1],
                  triplet[2][2] - triplet[0][2]]]

        # Finds vector perpendicular to the plane of the three points
        v = numpy.matrix([numpy.linalg.det([[trans[1][1], trans[2][1]],
                                            [trans[1][2], trans[2][2]]]),
                          numpy.linalg.det([[trans[1][0], trans[2][0]],
                                            [trans[1][2], trans[2][2]]]),
                          numpy.linalg.det([[trans[1][0], trans[2][0]],
                                            [trans[1][1], trans[2][1]]])])

        # Sets up first rotation matrix about the x axis
        theta = math.atan2(v.item(1), v.item(2))
        xmatr = numpy.matrix([[1, 0, 0], [0, math.cos(theta), - math.sin(theta)],
                              [0, math.sin(theta), math.cos(theta)]])
        trans1 = numpy.matrix(trans)
        rot1 = numpy.matrix.dot(trans1, xmatr)
        v1 = numpy.matrix.dot(v, xmatr)

        # Second rotation matrix about the y axis
        rho = math.atan2(v1.item(0), v1.item(2))
        ymatr = numpy.matrix([[math.cos(rho), 0, math.sin(rho)], [0, 1, 0],
                              [-math.sin(rho), 0, math.cos(rho)]])
        rot2 = numpy.matrix.dot(rot1, ymatr)

        # Should be in the xy plane now. Have to rotate such that two points
        # are on the x axis
        alph = math.atan2(rot2.item(4), rot2.item(3))
        bet = math.atan2(rot2.item(7), rot2.item(6))
        r1 = math.sqrt(math.pow(rot2.item(3), 2) + math.pow(rot2.item(4), 2))
        r2 = math.sqrt(math.pow(rot2.item(6), 2) + math.pow(rot2.item(7), 2))
        x = r1 / 2
        y = r2 * (1 - math.cos(bet - alph)) / (2.0 * math.sin(bet - alph))
        z = math.sqrt(abs(math.pow(dist, 2) - math.pow(x, 2) - math.pow(y, 2)))
        si_pos = numpy.matrix([x, y, z])

        # Rotate back to originial position
        init = math.atan2(si_pos.item(1), si_pos.item(0))
        r = math.sqrt(math.pow(si_pos.item(0), 2) + math.pow(si_pos.item(1), 2))
        x = r * math.cos(init + alph)
        y = r * math.sin(init + alph)
        si_pos = numpy.matrix([x, y, z])

        # Undo second rotation matrix
        iymatr = numpy.linalg.inv(ymatr)
        si_pos = numpy.matrix.dot(si_pos, iymatr)

        # Undo first rotation matrix
        ixmatr = numpy.linalg.inv(xmatr)
        si_pos = numpy.matrix.dot(si_pos, ixmatr)

        # Translate back so there is no point at the origin
        si_pos = [si_pos.item(0) + triplet[0][0],
                  si_pos.item(1) + triplet[0][1],
                  si_pos.item(2) + triplet[0][2]]

        return si_pos

        
            
class Center():
    """ An object that represents each individual center and keeps track of its
        information, such as its coordinates, ring size, and the distance
        between it and the hole """

    def __init__(self, coord, image_class):
        """ Constructor """
        
        self._coord = coord  # Coordinate where the center exists in pixels
        self._image = image_class
        
        self._nm_coord = self._image.pixelsToNm(self._coord)
        
        self._ring_size = ''  # Numbers ranging from 4 to 9
        self._hole_dist = ''  # In nm
        self._ring_center_coord = ''  # In pixels
        self._ring_center_nm_coord = ''  # In nm
        
    
    def getCoords(self):
        return self._coord
    
    
    def getRingSize(self):
        return self._ring_size
    
    
    def getHoleDist(self):
        return self._hole_dist
    
    
    def getRingCenterCoord(self):
        return self._ring_center_coord
    
    
    def getRingCenterNmCoord(self):
        return self._ring_center_nm_coord
    
    
    def assignRingSize(self, ring_size):
        self._ring_size = ring_size
    
    
    def assignHoleDist(self, dist):
        self._hole_dist = dist
        
    
    def assignRingCenterCoord(self, coord):
        self._ring_center_coord = coord
        self._ring_center_nm_coord = self._image.pixelsToNm(coord)



class Oxygen():
    """ A class for each oxygen atom found with its location and if it is a
        posible triplet """
    
    def __init__(self, coord):
        """ Constructor """
        
        self._coord = coord
        
        self._triplet = []  # A list of other oxygens objects close enough to
                            # it to have a single Si atom between them
        self._triplet_positions = []  # A list of the best triplet positions
        
    
    def getCoords(self):
        return self._coord
    
    
    def getTriplet(self):
        return self._triplet
    
    
    def getTripletPositions(self):
        return self._triplet_positions
    
    
    def addTriplet(self, atom):
        self._triplet.append(atom)
    
    
    def setTripletPositions(self, triplet_positions):
        self._triplet_positions = triplet_positions



def main():
    """ Initializes necessary objects """
    
    UserInterface() # Create the starting window and process initial information
    


if __name__ == "__main__":
    main()










