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


class STM():
    """ A class to describe the STM image. Includes information like filename,
    Image Dimensions (pixels), sample dimensions (nm), scale, number of holes,
    and coordinates of those holes."""
    
    def __init__(self):
        """ Constructor - gather preliminary info inputted into user interface """
        
        self._filename = Tk.StringVar()
        self._width = Tk.StringVar()
        self._height = Tk.StringVar()
        self._num_holes = Tk.StringVar()
        self._xyz_file = Tk.StringVar()
        
        # Image filename
        self._prompt = Tk.Label(master=root, text='Enter the filename of the image')
        self._prompt.pack(side=Tk.TOP)
        self._filenameEntry = Tk.Entry(master=root, textvariable=self._filename)
        self._filenameEntry.pack(side=Tk.TOP)
    
        # Width in nm
        self._widthPrompt = Tk.Label(master=root, text='Image Width (nm)')
        self._widthPrompt.pack(side=Tk.TOP)
        self._widthEntry = Tk.Entry(master=root, textvariable=self._width)
        self._widthEntry.pack(side=Tk.TOP)
    
        # Height in nm
        self._heightPrompt = Tk.Label(master=root, text='Image Height (nm)')
        self._heightPrompt.pack(side=Tk.TOP)
        self._heightEntry = Tk.Entry(master=root, textvariable=self._height)
        self._heightEntry.pack(side=Tk.TOP)
    
        # Number of holes
        self._holesPrompt = Tk.Label(master=root, text='Number of Holes')
        self._holesPrompt.pack(side=Tk.TOP)
        self._holesEntry = Tk.Entry(master=root, textvariable=self._num_holes)
        self._holesEntry.pack(side=Tk.TOP)
    
        # XYZ file
        self._xyzPrompt = Tk.Label(master=root, text='Enter XYZ Filename')
        self._xyzEntry = Tk.Entry(master=root, textvariable=self._xyz_file)

        # XYZ checkbox button
        self._xyzImport = Tk.IntVar()
        self._importXyzBtn = Tk.Checkbutton(master=root, text='Import XYZ ',
                                      variable=self._xyzImport,
                                      command=self.importXyz())
        self._importXyzBtn.pack(side=Tk.TOP)
    
        # Continue button
        self._continueBtn = Tk.Button(master=root, text='Continue', command=self.enterFilename())
        self._continueBtn.pack(side=Tk.TOP)
    
        # Keeps the window open and listens for events
        while True:
            try:
                root.mainloop()
                break
            # Fixes issue where scrolling accidentally would crash Tkinter
            except UnicodeDecodeError:
                pass
        
        # self._im_dim = im_dim  # [image width, image height] (pixels)
        # self._scale = im_dim[0] / sample_dim[0]  # ratio pixels/nm
        # self._hole_coords = []
        # self._sample_dim = initVars[1] # [sample width, sample height] (nm)


    def importXyz(self):
        if self._xyz_file:
            self._xyzPrompt.pack(side=Tk.TOP)
            self._xyzEntry.pack(side=Tk.TOP)
        else:
            self._xyzPrompt.pack_forget()
            self._xyzEntry.pack_forget()


    def enterFilename(self):
            # Get Filename from Text Entry
    
            try:
                io.imread(self._filename)
    
                # Destroy the entry UI
                self._prompt.destroy()
                self._filenameEntry.destroy()
                self._continueBtn.destroy()
                self._widthPrompt.destroy()
                self._heightPrompt.destroy()
                self._widthEntry.destroy()
                self._heightEntry.destroy()
                self._holesPrompt.destroy()
                self._holesEntry.destroy()
                self._xyzPrompt.destroy()
                self._xyzEntry.destroy()
                self._importXyzBtn.destroy()
    
                # Call the rest of the program using given parameters
                # centerFinder(filename, dimensions, num_holes, import_xyz,
                             xyz_filename)
            except ValueError:
                messagebox.showerror("Error", "Invalid dimensions and/or number of holes")
            except FileNotFoundError:
                messagebox.showerror("Error", "A file could not be found with that filename")


def main():
    """ Initializes all necessary objects """
    
    STM()
    


if __name__ == "__main__":
    main()










