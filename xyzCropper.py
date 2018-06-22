"""xyzCropper
   Will chop a certain amount off the left side of an image.
   call python3 xyzChopper.py filename cropDistance newfilename. """
import sys


def getCentersFromXyz(xyz_filename):
    centers = []
    with open(xyz_filename, encoding='utf-8') as f:
        file_lines = f.readlines()

    file_lines = [x.strip() for x in file_lines]
    for line in file_lines:
        split_line = line.split(" ")
        type_and_coord = [int(split_line[0]), float(split_line[1]) * 10, float(split_line[2]) * 10]
        centers.append(type_and_coord)
    return centers


def xyzCropper(filename, cropDistance):
    """crops xyz file"""
    centers = getCentersFromXyz(filename)
    goodCenters = []
    for center in centers:
        if center[1] > float(cropDistance):
            goodCenters.append(center)
    return goodCenters


def createXyzFile(centers, newFilename):
    text_file = open(newFilename[:-4], "w")
    for center in centers:
        ring_type = str(center[0])
        x_coord = str(center[1] / 10)
        y_coord = str(center[2] / 10)
        text_file.write(ring_type + ' ' + x_coord + ' ' + y_coord + ' 0\n')

    text_file.close()


def main():
    centers = xyzCropper(sys.argv[1], sys.argv[2])
    createXyzFile(centers, sys.argv[3])


main()
