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


class Startup():
    """ A starting class that gathers the preliminary information inputted in
        the user interface """
    
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
                                      command=self.importXyz())
        self._import_xyz_btn.pack(side=Tk.TOP)
    
        # Continue button
        self._continue_btn = Tk.Button(master=root, text='Continue', command=self.enterFilename())
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
        if self._xyz_file:
            self._xyz_prompt.pack(side=Tk.TOP)
            self._xyz_entry.pack(side=Tk.TOP)
        else:
            self._xyz_prompt.pack_forget()
            self._xyz_entry.pack_forget()


    def enterFilename(self):
            """ Get Filename from text entry """
    
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
                Image(self._filename, self._width, self._height,
                      self._num_holes, self._xyz_import, self._xyz_file)

            except ValueError:
                messagebox.showerror("Error", "Invalid dimensions and/or number of holes")
            except FileNotFoundError:
                messagebox.showerror("Error", "A file could not be found with that filename")



class Image():
    """ An overarching class to describe the STM image. Includes information
        like filename, image dimensions (pixels), sample dimensions (nm),
        scale, number of holes, and coordinates of those holes."""

    def __init__(self, filename, width, height, num_holes, xyz_import, xyz_file):
        """ Constructor """

        self._filename = filename
        self._sample_dim = [width, height]  # (nm)
        self._num_holes = num_holes
        self._xyz_import = xyz_import
        self._xyz_file = xyz_file
        
        self._image = self.importAndScale()
        
        self.grey()
        
        self._im_dim = (len(self._grey), len(self._grey[0])) # (image height, image width) (pixels)
        self._scale = self._im_dim[1] / self._sample_dim[0]  # ratio pixels/nm

        self._holes = Holes(self._num_holes, self._grey, self._im_dim, self)

        self._rings = Rings(self._num_holes, self._xyz_import, self._xyz_file,
                            self._grey_inv, self)
            

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
    
    
    def nmToPixels(self, nm_coord):
        """ Converts a set of nm coordinates into pixel coordinates """
        return [int(nm_coord[0] * self._scale), int(nm_coord[1] * self._scale)]
    
    
    def pixelsToNm(self, pixel_coord):
        """ Converts a set of pixel coordinates into nm coordinates """
        return [pixel_coord[0] / self._scale, pixel_coord[1] / self._scale]
    
    
    def pixelDistToNm(self, dist):
        """Converts a distance in pixels into a distance in nm """
        return dist / self._scale
    
    
    def returnHoleImage(self):
        """ Returns the mask image of the hole """
        
        return self._hole.returnHoleImage()
    
    
    def getNearestNeighbors(self, base_coords, search_coords, num_neighbors):
        """ Finds the distances and indices of the nearest given number of base
            coords to each of the provided search coords. """

        nearest = NearestNeighbors(n_neighbors=num_neighbors,
                                   algorithm='ball_tree').fit(base_coords)
        dist, ind = nearest.kneighbors(search_coords)

        return dist, ind
    
    
    def plotCirclesOnImage(self, image, coords, radius, color, reverse=False):
        """ Plots circles of a given radius and color on a given image at given coordinates """
        if not reverse:
            for coord in coords:
                rr, cc = draw.circle(coord[1], coord[0], radius, shape=self._im_dim)
                image[rr, cc] = color
        else:
            # Same as above but coord reversed
            for coord in coords:
                rr, cc = draw.circle(coord[0], coord[1], radius, shape=self._im_dim)
                image[rr, cc] = color
    
    

class Holes():
    """ A class to keep track of hole related information """
    
    def __init__(self, num_holes, grey, image_dim, image_class):
        """ Constructor """
        
        self._num_holes = num_holes
        self._grey = grey  # Grey scale image
        self._im_dim = image_dim # (image height, image width) (pixels)
        self._im_class = image_class
        
        self._hole_coords = self.getHoleCoords()
        self._hole_image = self.getHoleImage()  # Mask image of holes
    
    
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
    
    
    def returnHoleImage(self):
        """ Returns the mask image of the hole """
        
        return self._hole_image



class Rings():
    """ A class delaying with interpreting the rings within the STM image """

    def __init__(self, num_holes, xyz_import, xyz_file, grey_inv, image_class):
        """ Constructor """
        
        self._num_holes = num_holes
        self._xyz_import = xyz_import
        self._xyz_file = xyz_file
        self._grey_inv = grey_inv
        self._image = image_class
        
        # If not importing an XYZ file, identify Si centers on greyscale image
        if not self._xyz_import:
            self._imageblackOutHoles(self._image.returnHoleImage())
            self._centers = self.getCentersFromImage()
        else:
            self._centers = self.getCentersFromXYZ()
    
    
    def blackOutHoles(self, holes):
        """ Black out the holes in the greyscale image so no centers will be
            placed there """

        self._grey_inv = self._grey_inv * (1 - holes)
        self._grey_inv_copy = self._grey_inv

        self._peak_coord = feature.peak_local_max(self._grey_inv, min_distance=24,threshold_rel=0.2)
        
    
    def getCentersFromImage(self):
        """ Returns a list of the centers of the rings from the greyscale image"""
        
        blobs = feature.blob_dog(self._grey_inv, min_sigma=0.07, max_sigma=15,
                                 sigma_ratio=2.8, threshold=0.5, overlap=0.3)

        centers = self.getBlobCenters(blobs)
        
        # Find the average distance to closest neighbor
        c_dist, c_ind = self._image.getNearestNeighbors(centers, centers, 2)
        avg_closest = numpy.median(c_dist[:][1])
        num_iter = 2
        
        for i in range(num_iter):
            #Blackout regions around already found ring centers
            #TO DO: Possibly look into scaling this based on average closest
            average_thresh = 1.6

            self._image.plotCirclesOnImage(self._grey_inv, centers, avg_closest*average_thresh, 0)

            #Find blobs that may be rings that have not been found on first pass
            blobs = numpy.concatenate((blobs, feature.blob_dog(self._grey_inv, min_sigma=0.08,
                                                              max_sigma=20, sigma_ratio=2.8,
                                                              threshold=0.8, overlap=0.3)))
            centers = self.getBlobCenters(blobs)

            self._image.plotCirclesOnImage(self._grey_inv_copy, self._peak_coord,
                                           avg_closest * average_thresh, 0, True)

            self._peak_coord = numpy.concatenate((self._peak_coord,
                                                  feature.peak_local_max(self._grey_inv_copy, 
                                                     min_distance=24, threshold_rel=0.2)))
                
        return centers
    
    
    def getBlobCenters(self, blobs):
        """ Returns a list of the centers of the blobs """
        
        # Find blobs in image in order to locate ring centers:
        # min_sigma -> increase to find more small rings
        # max_sigma -> increase to find more large rings
        # sigma_ratio -> used in calculation of rings
        # threshold -> decrease to detect less intense rings
        # overlap -> fraction of the blobs that are allowed to overlap with each other

        centers = []
        for blob in blobs:
            x_coord = int(blob[1])
            y_coord = int(blob[0])
            centers.append([x_coord, y_coord])

        return centers
    
    
    def getCentersFromXyz(self):
        """ Read the XYZ file to find the centers noted within """
        
        centers = []
        with open(self._xyz_file, encoding='utf-8') as f:
            file_lines = f.readlines()

        file_lines = [x.strip() for x in file_lines]
        for line in file_lines:
            split_line = line.split(" ")
            nm_coord = [float(split_line[1]), float(split_line[2])]
            pixel_coord = self._image.nmToPixels(nm_coord)
            centers.append(pixel_coord)

        return centers



def main():
    """ Initializes necessary objects """
    
    Startup() # Create the starting window and process initial information
    


if __name__ == "__main__":
    main()










