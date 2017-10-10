
import math
import numpy

from skimage import io
from skimage import feature
from skimage import draw
from skimage import util
from skimage import color
from skimage import viewer
from skimage import morphology
from skimage import filters
from skimage import measure
from skimage import transform
from skimage import exposure

from sklearn.neighbors import NearestNeighbors

from scipy import ndimage as ndi

import matplotlib.pyplot as plt


filename = "170908_SiO2@RuVII_STM110_det (27.4x11.5).jpg"

light_centers = False

#Dimensions of image [width, height] in nm
dimensions = [27.4, 11.5]
#Number of holes in the image
num_holes = 1

#Outputs
save_image = False
output_filename = ""
save_xyz_file = False
xyz_output_filename = "Output.txt"

manual_correct = False
display_result = False

#Colors for each of the ring numbers
colors = [[178, 112, 248], [75, 176, 246], [67, 196, 127], [249, 222, 62], [249, 76, 62], [247, 38, 232]]

scaling_factor = 1

#View image in separate viewer
def viewImage(image):
    view = viewer.ImageViewer(image)
    view.show()
    return view

#Imports and scales the image by a given scaling factor
def importAndScale(filename):
    original = io.imread(filename)
    return transform.rescale(original, scaling_factor, mode='constant', preserve_range=True).astype('uint8')

image = importAndScale(filename)
image_copy = importAndScale(filename)

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

#Get coordinates of borders of holes from greyscale image
def getHoleCoords(greyscale, num_holes):
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

#Takes the coordinates of the borders of the holes and fills them in to get
# a mask image of the holes
def getHoleImage(hole_coords):
    hole_image = numpy.zeros((image_width,image_height))
    
    for coord in hole_coords:
        hole_image[coord[1],coord[0]] = 1

    return ndi.binary_fill_holes(morphology.binary_closing(hole_image))

hole_coords = getHoleCoords(grey, num_holes)

holes = getHoleImage(hole_coords)

#viewImage(holes)

#Black out the holes in the greyscale image so no centers will be placed there
grey_inv = grey_inv * (1 - holes)

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

#blobs = feature.blob_dog(grey_inv, min_sigma=0.07, max_sigma=30, sigma_ratio=2.8, threshold=0.55, overlap=0.26)
blobs = feature.blob_dog(grey_inv, min_sigma=0.07, max_sigma=30, sigma_ratio=2.8, threshold=0.59, overlap=0.3)
centers = getBlobCenters(blobs)

#Find the average distance to closest neighbor
c_dist, c_ind = getNearestNeighbors(centers, centers, 2)
average_closest = numpy.median(c_dist[:][1])

#Blackout regions around already found ring centers
average_thresh = 1.4

plotCirclesOnImage(grey_inv, centers, average_closest*average_thresh, 0)

#Find blobs that may be rings that have not been found on first pass
blobs = numpy.concatenate((blobs,feature.blob_dog(grey_inv, min_sigma=1, max_sigma=30, sigma_ratio=2.8, threshold=1.1, overlap=0.5)))
centers = getBlobCenters(blobs)

#Converts pixel coordinates to nm coordinates based on the dimensions of the image
def pixelsToNm(pixel_coord, nm_dim, im_width, im_height):
    x_scale = nm_dim[0]/im_width
    y_scale = nm_dim[1]/im_height
    
    return [pixel_coord[0]*x_scale, pixel_coord[1]*y_scale]

def onclick(event):
    #Remove the center that was clicked
    if numpy.array_equal(image_copy[int(event.ydata),int(event.xdata)], [0, 0, 250]):
        min_dist = math.hypot(event.xdata-centers[0][0], event.ydata-centers[0][1])
        match_ind = 0
        for i in range(len(centers)):
            cur_dist = math.hypot(event.xdata-centers[i][0], event.ydata-centers[i][1])
            if cur_dist < min_dist:
                min_dist = cur_dist
                match_ind = i
        del centers[match_ind]
    #Add a center where the user clicked
    else:
        centers.append([event.xdata, event.ydata]);

def manualCenterCorrect(image_copy, centers):
    #Plot preliminary center locations on image
    plotCirclesOnImage(image_copy, centers, 10, [0, 0, 250])

    point_select = viewer.ImageViewer(image_copy)
    point_select.canvas.mpl_connect('button_press_event', onclick)

    point_select.show()
    
if manual_correct:
    manualCenterCorrect(image_copy, centers)

def getNumNeighbors(centers, thresh, average_closest):
    #Get distances and indices of 9 nearest neighbors to every center
    distances, indices = getNearestNeighbors(centers, centers, 10)

    #Gets the distances from every center to the nearest point on the edge of a hole
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
    
        exclude_thresh = 1.6
    
        #Gets the coordinates of the circles around the centers
        if 4 <= scaled_num_neighbors <= 9 and percent_visible < 1.8 and hole_distances[k][1] > exclude_thresh*average_closest:
            hole_dist.append(hole_distances[k][1])
            ring_size.append(scaled_num_neighbors)
            center_coord.append(centers[k])
        
    return hole_dist, ring_size, center_coord

#Threshold for the maximum distance that two centers can be apart to be concidered neighbors
thresh = 1.35
    
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
    text_file = open(xyz_output_filename, "w")
    
    for i in range(len(ring_size)):
        text_file.write(str(ring_size[i]) + ' '+ str(center_coord[i][0]) + ' ' + str(center_coord[i][1]) + ' 0\n')
    
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

def plotBinHist(bin_list):
    #Histograms are normalized so they can be compared on same scale
    plt.hist(bin_list, bins='auto', normed=True)
    
def plotRingSizePercent(bin_list, bin_mids):    
    size_perc = [[],[],[],[],[],[]]
    
    for cur_bin in bin_list:
        bin_perc = [0,0,0,0,0,0]
        
        for size in cur_bin:
            bin_perc[size-4] += 1
            
        bin_perc = [x / len(cur_bin) for x in bin_perc]
        
        for k in range(6):
            size_perc[k].append(bin_perc[k])
        
    for i in range(len(size_perc)): 
        plt.plot(bin_mids, size_perc[i], label=str(i+4) + ' MR')
        
    legend = plt.legend(loc=0, ncol=2)
    plt.gca().add_artist(legend)
    
    plt.show()

if save_xyz_file:
    createXyzFile(center_coord, ring_size)

bin_list, bin_mids = splitRingsIntoBins(180, hole_dist, ring_size)

#plotBinHist(bin_list)

plotRingSizePercent(bin_list, bin_mids)

if display_result:
    #Display image with centers plotted
    viewImage(image)
    
    if save_image:
        io.imsave(output_filename, image)
