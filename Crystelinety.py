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

light_centers = False

#170908_SiO2@RuVII_STM110_det (27.4x11.5).jpg

save_figures = False

#Colors for each of the ring numbers
colors = [[178, 112, 248], [75, 176, 246], [67, 196, 127], [249, 222, 62], [249, 76, 62], [247, 38, 232]]

scaling_factor = 1

#Sets up the Tk window
root = Tk.Tk()
root.wm_title("Auto Ring Center Finder")


bin_list = [[5,6,7,8],[6,6,7],[6,6,6,6,6,8],[6,6,6,6,6,6]]
bin_mids = [0.5,1,1.5,2]
def Crystallinity(bin_list, bin_mids):
    
    crys_list = []
    for cur_list in bin_list:
        num_sixes = cur_list.count(6)
        crystallinity = num_sixes/len(cur_list)
        crys_list.append(crystallinity)
                    
    plt.plot(bin_mids, crys_list, 'ro')
    plt.axis([0, bin_mids[-1] + 1, 0, 1])
    plt.title('Crystallinity')
    plt.xlabel('Distance')
    plt.ylabel('Crystallinity')
    plt.show()

    
    
Crystallinity(bin_list, bin_mids)