# -*- coding: utf-8 -*-
"""
Created on Tue May 30 09:40:34 2017

@author: Kristen
"""

import math
import numpy

def distance(position1,position2):
    return math.sqrt(math.pow(position1[0]-position2[0],2)
                     +math.pow(position1[1]-position2[1],2)
                     +math.pow(position1[2]-position2[2],2))

def locate_si(positions,dist):
    #assumes presorted positions by x position
    doubles=[]

    #finds all within the given radius and adds those doubles to the list
    for i in range(len(positions)):
        for j in range(1,len(positions)-i):
            if distance(positions[i],positions[i+j])<=dist:
                doubles.append([positions[i],positions[i+j]])

    return doubles

def find_o(positions,dist):
    
    opositions=[]
    
    for i in range(len(positions)):
        #center at origin
        pos1=[0,0,0]
        pos2=[positions[i][1][0]-positions[i][0][0],positions[i][1][1]-
              positions[i][0][1],positions[i][1][2]-positions[i][0][2]]


        #rotate until both points are in the xy plane
        theta=numpy.arctan2(pos2[1],pos2[0])
        phi=numpy.arctan2(pos2[2],pos2[0])
        newx=math.sqrt(math.pow(pos2[0],2)+math.pow(pos2[2],2))
        newy=newx*math.tan(theta)

        #find in si position (midpoint between origin and pos 2 in the x-y 
        # plane with x making up the difference)
        x=newx/2
        y=newy/2
        
        if math.pow(dist,2)-math.pow(x,2)-math.pow(y,2)>0:
            z=math.sqrt(math.pow(dist,2)-math.pow(x,2)-math.pow(y,2))
        else:
            z=0

        #current angle above x-y plane
        r=math.sqrt(math.pow(x,2)+math.pow(y,2)+math.pow(z,2))
        alph=math.asin(z/r)

        #when rotated back,it will rotate to angle phi+alph
        opos=[r*math.cos(theta)*math.cos(alph+phi),
              r*math.sin(theta)*math.cos(alph+phi),
                        r*math.sin(alph+phi)]

        #append to the list
        opositions.append([opos[0]+positions[i][0][0],opos[1]+
                           positions[i][0][1],opos[2]+positions[i][0][2]])

    return opositions

def main():

    #convert data in file into floats and append to a position list
    with open("SiPositionsExample.txt") as f:
        content=f.readline()
    
    string=""

    locations=[]

    for i in range(len(content)):
        if content[i] == " ":
            locations.append(float(string))
            string=""
        else:
            string+=content[i]

    locations.append(float(string))

    positions=[[""]]

    for i in range(len(locations)):
        if i%3==0:
            positions[i/3]=[locations[i]]
            positions.append("")
        else:
            positions[i/3].append(locations[i])

    del positions[len(positions)-1]

    #sort positions for the double finder function

    positions=sorted(positions)

    #find doubles
    doubles=locate_si(positions,.36)

    remov=[]

    for i in range(len(doubles)):
        if len(doubles[i])!=2:
            remov.insert(0,i)

    for i in range(len(remov)):
        del doubles[remov[i]]

    #find O positions
    o_locations=find_o(doubles,.17)

    #write Si positions to an out file
    out=open("O Positions Output_2.txt","w")
    out.write(str(o_locations))
    out.write("nn")

if __name__== "__main__":
    main()

