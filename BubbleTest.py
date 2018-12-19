import math
import numpy
import matplotlib

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


filename = 'AmorpHole.jpg'

image = io.imread(filename)
grey = color.rgb2grey(image)
grey = exposure.equalize_adapthist(exposure.adjust_gamma(grey))
# grey = exposure.equalize_adapthist(exposure.adjust_gamma(grey))
grey = grey > 0.5

# grey = util.invert(grey)


# def getHoleImage(hole_coords):
#         hole_image = numpy.zeros((image_height,image_width))
        
#         for coord in hole_coords:
#             hole_image[coord[1],coord[0]] = 1

#         return ndi.binary_fill_holes(morphology.binary_closing(hole_image))


# hole_coords = getHoleCoords(grey, 1)
# holes = getHoleImage(hole_coords)



io.imshow(grey)
io.show()