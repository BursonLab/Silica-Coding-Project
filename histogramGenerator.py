"""
Elisabeth's pair distance and g(r) histogram generator

Primarily for qualitative comparison of amorphous and crystaline bubble rafts
g(r) defined as pair distance divided by the area covered by the ring
Edge effects entirely ignored

4/14/19 - 5/11/19

References:
1. http://www.physics.emory.edu/faculty/weeks//idl/gofr.html
2. https://jakevdp.github.io/PythonDataScienceHandbook/04.05-
    histograms-and-binnings.html
"""

import numpy as np
import math
from matplotlib import pyplot as plt


def radius_calculator(centers):
    """Given a list of centers, makes a list of the distances between the
       centers."""
    distances = []

    # Iterates through each center in comparison to all other centers
    for a in centers:
        for b in centers:
            # Ensures each radius is only counted once
            if centers.index(a) < centers.index(b):
                radius = ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** .5
                distances.append(radius)

    return distances


def pair_distance_histogram(data):
    """Given a list of distances, creates a historgram."""

    # Takes user input on how many bins to create
    bin_num = (int(input("Enter the desired number of bins in the histogram:"
                         )) + 1)

    # Determines the numeric locations of the bins
    bins = np.linspace(math.floor(min(data)), math.ceil(max(data)), bin_num)

    # Plots the histogram
    plt.hist(data, bins=bins)
    plt.title('Pair Distance Histogram')
    plt.xlabel('Radius')
    plt.ylabel('Frequency')

    plt.show()


def g_r_histogram(data):
    """Given a list of distances, creates a historgram divided by radius"""

    # Takes user input on how many bins to create
    bin_num = (int(input("Enter the desired number of bins in the histogram:"
                         )) + 1)

    # Calculates the histogram without plotting it
    counts, bins = np.histogram(data, bins=bin_num)

    bi = list(bins)
    bi.pop()

    # Updates the counts to account for area
    new_counts = [counts[i] / (math.pi * (bins[i + 1] ** 2 - bins[i] ** 2))
                  for i in range(len(counts))]

    # Creates a bar plot of the counts per bin (essentially a histogram)
    plt.bar(bi, new_counts, width=bi[1] - bi[0], align='edge')
    plt.title('g(r) Histogram')
    plt.xlabel('Radius (pixels)')
    plt.ylabel('Frequency / Area')
    plt.show()


def main():

    # Gets the list of centers from a file
    file_name = input("Enter the name of the file with the list of centers: ")
    file = open(file_name)
    centers = eval(file.read())

    # Calculates the radii
    radii = radius_calculator(centers)

    # Asks user to choose the desired histogram type
    which_hist = input("Enter the type of histogram you wish to"
                       " see (pair distance or g(r)): ")

    if which_hist == "pair distance":
        pair_distance_histogram(radii)

    elif which_hist == "g(r)":
        g_r_histogram(radii)


if __name__ == '__main__':
    main()
