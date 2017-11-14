
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

light_centers = False

#170908_SiO2@RuVII_STM110_det (27.4x11.5).jpg

save_figures = False

#Colors for each of the ring numbers
colors = [[178, 112, 248], [75, 176, 246], [67, 196, 127], [249, 222, 62], [249, 76, 62], [247, 38, 232]]

scaling_factor = 1

#Sets up the Tk window
root = Tk.Tk()
root.wm_title("Auto Ring Center Finder")

def getFilename():
    
    #Variables for storing text entered
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
        #Get Filename from Text Entry
        filename = entryFilename.get()
        xyz_filename = entryXyz.get()
        import_xyz = xyzImport.get()
        
        try:
            io.imread(filename)
            
            dimensions = [float(entryWidth.get()), float(entryHeight.get())]
            
            num_holes = int(entryHoles.get())
            
            #Destroy the entry UI
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
        
            #Call the rest of the program using given parameters
            centerFinder(filename, dimensions, num_holes, import_xyz, xyz_filename)
        except ValueError:
            messagebox.showerror("Error", "Invalid dimensions and/or number of holes")
        except FileNotFoundError:
            messagebox.showerror("Error", "A file could not be found with that filename")
        
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
    
    importXyzBtn = Tk.Checkbutton(master=root, text='Import XYZ ', variable=xyzImport, command=importXyz)
    importXyzBtn.pack(side=Tk.TOP)
    
    continueBtn = Tk.Button(master=root, text='Continue', command=enterFilename)
    continueBtn.pack(side=Tk.TOP)
    
    
    
    #Keeps the window open and listens for events
    while True:
        try:
            root.mainloop()
            break
        #Fixes issue where scrolling accidentally would crash Tkinter
        except UnicodeDecodeError:
            pass
        

def centerFinder(filename, dimensions, num_holes, import_xyz, xyz_filename):   
        
    #Converts pixel coordinates to nm coordinates based on the dimensions of the image
    def pixelsToNm(pixel_coord, nm_dim, im_width, im_height):
        x_scale = nm_dim[0]/im_width
        y_scale = nm_dim[1]/im_height
        
        return [pixel_coord[0]*x_scale, pixel_coord[1]*y_scale]
    
    def nmToPixels(nm_coord, nm_dim, im_width, im_height):
        x_scale = im_width/nm_dim[0]
        y_scale = im_height/nm_dim[1]
        
        return [int(nm_coord[0]*x_scale), int(nm_coord[1]*y_scale)]
    
        
    #Imports and scales the image by a given scaling factor
    def importAndScale(filename):
        original = io.imread(filename)
        return transform.rescale(original, scaling_factor, mode='constant', preserve_range=True).astype('uint8')
    
    image = importAndScale(filename)
    
    #Convert image to greyscale 
    grey = color.rgb2grey(image)
    
    #Adjust the brightness of the greyscale to make image more uniform
    grey = exposure.equalize_adapthist(exposure.adjust_gamma(grey))
    
    #Invert image if centers are dark, else keep same
    if light_centers:
        grey_inv = grey
    else:
        grey_inv = util.invert(grey)
    
    image_width = len(grey)
    image_height = len(grey[0])
    
    #Plots circles of a given radius and color on a given image at given coordinates
    def plotCirclesOnImage(image, coords, radius, color):
        for coord in coords:
            rr, cc = draw.circle(coord[1], coord[0], radius, shape=(image_width, image_height))
            image[rr, cc] = color
            
            
    #Finds the distances and indices of the nearest given number of the base coords
    # to each of the provided search coords
    def getNearestNeighbors(base_coords, search_coords, num_neighbors):
        nearest = NearestNeighbors(n_neighbors=num_neighbors, algorithm='ball_tree').fit(base_coords)
        dist, ind = nearest.kneighbors(search_coords)
        return dist, ind
    

    #Get coordinates of borders of holes from greyscale image
    def getHoleCoords(greyscale, num_holes):
        if num_holes > 0:
            #Threshold greyscale image to create a binary image
            im_thresh = filters.threshold_minimum(greyscale)
            binary = greyscale > im_thresh
        
            #Erode the binary image and then subtract from original to get borders
            erosion = numpy.pad(morphology.binary_erosion(binary)[2:-2,2:-2],2,'maximum')
            borders = binary ^ erosion
            borders = numpy.pad(borders[1:-1,1:-1],1,'edge')
            
            #Label the regions in the image
            label_image = measure.label(borders)
            regions = measure.regionprops(label_image)
        
            #Put all of the areas of regions in border image into list
            areas = []
            for region in regions:
                areas.append(region.area)
    
            #Get the coordinates of the largest num_holes number of border regions
            hole_coords = []
            for i in range(num_holes):
                max_ind = numpy.argmax(areas)
                coords = regions[max_ind].coords
                for coord in coords:
                    hole_coords.append([coord[1],coord[0]])
                
                del regions[max_ind]
                del areas[max_ind]
        
            return hole_coords
        else:
            return []
    
    #Takes the coordinates of the borders of the holes and fills them in to get
    # a mask image of the holes
    def getHoleImage(hole_coords):
        hole_image = numpy.zeros((image_width,image_height))
        
        for coord in hole_coords:
            hole_image[coord[1],coord[0]] = 1
    
        return ndi.binary_fill_holes(morphology.binary_closing(hole_image))
    
    hole_coords = getHoleCoords(grey, num_holes)
    
    holes = getHoleImage(hole_coords)
    
    if not import_xyz:
        #Black out the holes in the greyscale image so no centers will be placed there
        grey_inv = grey_inv * (1 - holes)
        
        #Returns a list of the centers of the blobs
        def getBlobCenters(blobs):
            centers = []
            
            for blob in blobs:
                x_coord = int(blob[1])
                y_coord = int(blob[0])
            
                centers.append([x_coord, y_coord])
            
            return centers
        
        #Find blobs in image in order to locate ring centers
    
        # min_sigma -> increase to find more small rings
        # max_sigma -> increase to find more large rings
        # sigma_ratio -> used in calculation of rings
        # threshold -> decrease to detect less intense rings
        # overlap -> fraction of the blobs that are allowed to overlap with each other
    
        #blobs = feature.blob_dog(grey_inv, min_sigma=0.03, max_sigma=30, sigma_ratio=2.8, threshold=0.8, overlap=0.5)
        blobs = feature.blob_dog(grey_inv, min_sigma=0.07, max_sigma=15, sigma_ratio=2.8, threshold=0.57, overlap=0.3)
        centers = getBlobCenters(blobs)
        
        #Find the average distance to closest neighbor
        c_dist, c_ind = getNearestNeighbors(centers, centers, 2)
        avg_closest = numpy.median(c_dist[:][1])
        
        num_iter = 1
        
        for i in range(num_iter):
            #Blackout regions around already found ring centers
            #TODO: Possibly look into scaling this based on average closest 
            average_thresh = 1.6 #2.1#1.6 #1.92
        
            plotCirclesOnImage(grey_inv, centers, avg_closest*average_thresh, 0)
        
            #Find blobs that may be rings that have not been found on first pass
            blobs = numpy.concatenate((blobs,feature.blob_dog(grey_inv, min_sigma=0.08, max_sigma=20, sigma_ratio=2.8, threshold=0.8, overlap=0.3)))
            centers = getBlobCenters(blobs)
    
    else:
        
        def getCentersFromXyz(xyz_filename):
            centers = []
            
            with open(xyz_filename, encoding='utf-8') as f:
                file_lines = f.readlines()
    
            file_lines = [x.strip() for x in file_lines] 
            
            for line in file_lines:
                split_line = line.split(" ")
                
                nm_coord = [float(split_line[1]), float(split_line[2])]
                pixel_coord = nmToPixels(nm_coord, dimensions, image_width, image_height)
                
                centers.append(pixel_coord)
        
            return centers
        
        centers = getCentersFromXyz(xyz_filename)
        
        
    #Find the average distance to closest neighbor
    c_dist, c_ind = getNearestNeighbors(centers, centers, 2)
    average_closest = numpy.median(c_dist[:][1])
    
    def getNumNeighbors(centers, thresh, average_closest):
        #Get distances and indices of 9 nearest neighbors to every center
        distances, indices = getNearestNeighbors(centers, centers, 10)
    
        #Gets the distances from every center to the nearest point on the edge of a hole
        if num_holes > 0:
            hole_distances, hole_inds = getNearestNeighbors(hole_coords, centers, 2)
    
        hole_dist = []
        ring_size = []
        center_coord = []
    
        for k in range(len(distances)):
            n_dists = distances[k]
         
            #Averages the distances to closest 4 centers and multiplies by a 
            # threshold to get the max distance for something to be a neighbor
            max_dist = numpy.mean(n_dists[1:5])*thresh
            
            #Determines how many of the neighbors are within the max distance
            num_neighbors = 9
            for i in range(4,10):
                if n_dists[i] > max_dist:
                    num_neighbors = i-1
                    break
        
        
            r_full, c_full = draw.circle(centers[k][1], centers[k][0], max_dist);
            r_bound, c_bound = draw.circle(centers[k][1], centers[k][0], max_dist, shape=(image_width, image_height));
        
            #Gets the percentage of the ring neighbors that are visible in the window
            percent_visible = len(r_full)/len(r_bound)
        
            #Scales the number of neighbors based on what it should be if all the 
            # ring neighbors were visible
            scaled_num_neighbors = int(num_neighbors * percent_visible)
        
            exclude_thresh = 1.9
            
            if num_holes > 0:
                cur_hole_dist = hole_distances[k][1]
            else:
                cur_hole_dist = 2*exclude_thresh*average_closest
                
            #Gets the coordinates of the circles around the centers
            if 4 <= scaled_num_neighbors <= 9 and percent_visible < 1.2 and cur_hole_dist > exclude_thresh*average_closest:
                hole_dist.append(cur_hole_dist)
                ring_size.append(scaled_num_neighbors)
                center_coord.append(centers[k])
            
        return hole_dist, ring_size, center_coord
    
    #Threshold for the maximum distance that two centers can be apart to be concidered neighbors
    #thresh = 1.35
    thresh = 1.48
    #thresh = 1.52
        
    hole_dist, ring_size, center_coord = getNumNeighbors(centers, thresh, average_closest)
    
    def plotRingCenters(image, ring_size, centers, average_closest):
        for i in range(len(ring_size)):
            #Get circle coordinates for outlines and actual circles for centers
            r_out, c_out = draw.circle(centers[i][1], centers[i][0], int(average_closest/3)+3, shape=(image_width, image_height))
            rr, cc = draw.circle(centers[i][1], centers[i][0], int(average_closest/3), shape=(image_width, image_height))
                    
            #Plot outlines
            image[r_out, c_out] = [0, 0, 0]
            
            #Assign appropriate colors to center coordinates
            image[rr, cc] = colors[ring_size[i]-4]
            
    plotRingCenters(image, ring_size, center_coord, average_closest)
            
    #Outputs the center data in xyz file format
    def createXyzFile(center_coord, ring_size):
        text_file = open(filename[:-4]+'xyz.txt', "w")
        
        for i in range(len(ring_size)):
            nmCoord = pixelsToNm(center_coord[i], dimensions, image_width, image_height)
            text_file.write(str(ring_size[i]) + ' '+ str(nmCoord[0]) + ' ' + str(nmCoord[1]) + ' 0\n')
        
        text_file.close()
        
    #Splits the ring size and hole distance data into bins based on 
    # a given bin size and plots histograms of all of the bins 
    def splitRingsIntoBins(bin_size, hole_dist, ring_size):
        bin_list = []
        bin_mids = []
        
        max_dist = int(numpy.amax(hole_dist))
        
        for k in range(0, max_dist, bin_size):
            bin_mids.append(k + (bin_size/2))
            
            cur_bin = []
            for i in range(len(hole_dist)):
                if k <= hole_dist[i] < k + bin_size:
                    cur_bin.append(ring_size[i])
            
            bin_list.append(cur_bin)
        
        return bin_list, bin_mids
    """
    def plotBinHist(bin_list):
        #Histograms are normalized so they can be compared on same scale
        plt.hist(bin_list, bins='auto', normed=True)
    """
    
    def createWindow(image):
        print("Launched Auto Ring Finder")
        fig = Figure(figsize=(10, 6), dpi=100)
        ax = fig.add_subplot(111)
        ax.imshow(image)
        
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
    
        toolbar = NavigationToolbar2TkAgg(canvas, root)
        toolbar.update()
        canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        
        cid = []
        
        binSizeTxt = Tk.StringVar()
                        
        def addcenter(event):
            #Add a center where the user clicked
            centers.append([event.xdata, event.ydata]);
    
            #Replot image
            replotImage()
            
        def removecenter(event):
            #Remove the nearest center
            min_dist = math.hypot(event.xdata-centers[0][0], event.ydata-centers[0][1])
            match_ind = 0
            for i in range(len(centers)):
                cur_dist = math.hypot(event.xdata-centers[i][0], event.ydata-centers[i][1])
                if cur_dist < min_dist:
                    min_dist = cur_dist
                    match_ind = i
            del centers[match_ind]
            
            #Replot image
            replotImage()
            
        def movecenter(event):
            #Remove the nearest center
            min_dist = math.hypot(event.xdata-centers[0][0], event.ydata-centers[0][1])
            match_ind = 0
            for i in range(len(centers)):
                cur_dist = math.hypot(event.xdata-centers[i][0], event.ydata-centers[i][1])
                if cur_dist < min_dist:
                    min_dist = cur_dist
                    match_ind = i
            del centers[match_ind]
            
            #Replaces the center with one in the location clicked
            centers.append([event.xdata, event.ydata]);
            
            #Replot image
            replotImage()
            
        def replotImage():
            #Replot image
            image = io.imread(filename)
            
            hole_dist, ring_size, center_coord = getNumNeighbors(centers, thresh, average_closest)
            plotRingCenters(image, ring_size, center_coord, average_closest)
            
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
        
        addRingBtn = Tk.Button(master=root, text='Add Ring', command=ringAdd)
        addRingBtn.pack(side=Tk.LEFT)
        
        removeRingBtn = Tk.Button(master=root, text='Remove Ring', command=ringRemove)
        removeRingBtn.pack(side=Tk.LEFT)
        
        moveRingBtn = Tk.Button(master=root, text='Move Ring', command=ringMove)
        moveRingBtn.pack(side=Tk.LEFT)
        
        def saveImage():
            print("Saving image...")
            io.imsave(filename[:-4] + 'plotted.jpg', image)
        
        saveBtn = Tk.Button(master=root, text='Save Image', command=saveImage)
        saveBtn.pack(side=Tk.RIGHT)
        
        def plotRingSizePercent(bin_list, bin_mids): 
            ax = fig.add_subplot(111)
            
            size_perc = [[],[],[],[],[],[]]
        
            for cur_bin in bin_list:
                #Initializes the counts of the ring sizes for the bin at zero
                bin_perc = [0,0,0,0,0,0]
            
                #Increments the bin perc for each ring of that size
                for size in cur_bin:
                    bin_perc[size-4] += 1
            
                #Divides all ring size totals by the total number of rings in bin
                bin_perc = [x / len(cur_bin) for x in bin_perc]
            
                for k in range(6):
                    size_perc[k].append(bin_perc[k])
            
            #Plots the mid point of the bin vs. percentages of each ring size
            for i in range(len(size_perc)): 
                ax.plot(bin_mids, size_perc[i], label=str(i+4) + ' MR')
        
            #Adds a legend to the plot
            legend = ax.legend(loc=0, ncol=2)
            ax.add_artist(legend)
        
            #plt.ylabel('Percentage of Rings')
            #plt.xlabel('Distance from Nearest Hole (px)')
    
            canvas.draw()
            
        def percPlot():
            #See if the bin size is valid, and if so, plot the ring size percentages
            try:
                binSize = int(binSizeTxt.get())
                bin_list, bin_mids = splitRingsIntoBins(binSize, hole_dist, ring_size)
                fig.clf()
                plotRingSizePercent(bin_list, bin_mids)
            except ValueError:
                messagebox.showerror("Error", "Invalid Bin Size (must be an integer)")
                
        def xyzFile():
            hole_dist, ring_size, center_coord = getNumNeighbors(centers, thresh, average_closest)
            createXyzFile(center_coord, ring_size)
        
        def doneEditing():
            #Destroy editing buttons
            addRingBtn.destroy()
            removeRingBtn.destroy()
            moveRingBtn.destroy()
            doneBtn.destroy()
            
            fig.clf()
            canvas.draw()        
            
            #Label for bin size text entry
            binSizeLabel = Tk.Label(master=root, text='Bin Size')
            binSizeLabel.pack(side=Tk.LEFT)
            
            #Bin size text entry field
            binSizeEntry = Tk.Entry(master=root, textvariable=binSizeTxt)
            binSizeEntry.pack(side=Tk.LEFT)
            
            #Button to generate ring size percentage plot
            percPlotBtn = Tk.Button(master=root, text='Ring Size Percentage Plot', command=percPlot)
            percPlotBtn.pack(side=Tk.LEFT)
            
            #Button to export to an xyz file
            xyzFileBtn = Tk.Button(master=root, text='Create xyz File', command=xyzFile)
            xyzFileBtn.pack(side=Tk.LEFT)
            
        #Button to finish editing and move to exporting / plotting graphs    
        doneBtn = Tk.Button(master=root, text='Done Editing', command=doneEditing)
        doneBtn.pack(side=Tk.LEFT)
        
    
    createWindow(image)
    
getFilename()