
import math
import numpy

from skimage import io
from skimage import feature
from skimage import draw
from skimage import util
from skimage import color
from skimage import viewer

from sklearn.neighbors import NearestNeighbors


filename = "SampleImage.jpg"#"STM.jpg"#"connectionPoint.png"

#Dimensions of image [width, height] in nm
dimensions = [9.3, 5.6]

#Outputs
save_image = False
output_filename = "Output.png"
save_xyz_file = False
xyz_output_filename = "Output.txt"

#Key Image
has_key = False
key_filename = "connectionPointKey.png"

display_result = True

#Colors for each of the ring numbers
colors = [[178, 112, 248], [75, 176, 246], [67, 196, 127], [249, 222, 62], [249, 76, 62], [247, 38, 232]]

image = io.imread(filename)
image_copy = io.imread(filename)
key_image = io.imread(key_filename)

#Convert image to greyscale 
grey = color.rgb2grey(image)

#Invert image
grey_inv = util.invert(grey)

image_width = len(grey)
image_height = len(grey[0])

#Returns a list of the centers of the blobs
def getBlobCenters(blobs):
    centers = []
    
    for blob in blobs:
        x_coord = int(blob[1])
        y_coord = int(blob[0])
    
        centers.append([x_coord, y_coord])
    
    return centers
    
#Find blobs in image in order to locate ring centers
blobs = feature.blob_dog(grey_inv, min_sigma=0.03, max_sigma=30, sigma_ratio=2.8, threshold=0.8, overlap=0.5)

centers = getBlobCenters(blobs)

#Find the average distance to closest neighbor
closest_pair = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(centers)
c_dist, c_ind = closest_pair.kneighbors(centers)
average_closest = numpy.median(c_dist[:][1])

#Blackout regions around already found ring centers
average_thresh = 1.4
for center in centers:
    rr, cc = draw.circle(center[1], center[0], average_closest*average_thresh, shape=(image_width, image_height))
    grey_inv[rr, cc] = 0
    
#Find blobs that may be rings that have not been found on first pass
blobs = numpy.concatenate((blobs,feature.blob_dog(grey_inv, min_sigma=1, max_sigma=30, sigma_ratio=2.8, threshold=1.1, overlap=0.5)))

centers = getBlobCenters(blobs)

def pixelsToNm(pixel_coord, nm_dim, im_width, im_height):
    x_scale = nm_dim[0]/im_width
    y_scale = nm_dim[1]/im_height
    
    return [pixel_coord[0]*x_scale, pixel_coord[1]*y_scale]

for center in centers:
    rr, cc = draw.circle(center[1], center[0], 10, shape=(image_width, image_height))
    image_copy[rr, cc] = [0,0,250]


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
    
point_select = viewer.ImageViewer(image_copy)
cid = point_select.canvas.mpl_connect('button_press_event', onclick)

point_select.show()

#fig.show()

#Get distances and indices of 9 nearest neighbors to every center
neighbors = NearestNeighbors(n_neighbors=10, algorithm='ball_tree').fit(centers)
distances, indices = neighbors.kneighbors(centers)


thresh = 1.35
draw_lines = True

if save_xyz_file:
    text_file = open(xyz_output_filename, "w")

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
    
    #Gets the coordinates of the circles around the centers
    r_out, c_out = draw.circle(centers[k][1], centers[k][0], 12, shape=(image_width, image_height))
    rr, cc = draw.circle(centers[k][1], centers[k][0], 9, shape=(image_width, image_height))
    
    
    if scaled_num_neighbors <= 9 and percent_visible < 1.8:
        #Draws lines connecting a point to its neighbors
        if draw_lines:
            for i in range(1, scaled_num_neighbors):
                r_line, c_line = draw.line(centers[k][1], centers[k][0], centers[indices[k][i]][1], centers[indices[k][i]][0])
                try:
                    image[r_line, c_line] = [250, 250, 250]
                except IndexError:
                    pass
                
        image[r_out, c_out] = [0, 0, 0]
        
        #Assign appropriate colors to center coordinates
        image[rr, cc] = colors[scaled_num_neighbors-4]
        
        if save_xyz_file:
            text_file.write(str(scaled_num_neighbors) + ' '+ str(centers[k][0]) + ' ' + str(centers[k][1]) + ' 0\n')

if save_xyz_file:
    text_file.close()

if display_result:
    #Display image with centers plotted
    if has_key:
        image_view = viewer.CollectionViewer([image, key_image])
    else:
        image_view = viewer.ImageViewer(image)
    
    if save_image:
        io.imsave(output_filename, image)
    image_view.show()
