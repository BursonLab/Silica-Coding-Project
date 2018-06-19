# -*- coding: utf-8 -*-
"""
Created on Mon May 29 16:04:03 2017

@author: Kristen
"""
import math
import numpy


#finds the distance between two atoms
def distance(position1, position2):
    return math.sqrt(math.pow(position1[0] - position2[0],2) + 
                     math.pow(position1[1] - position2[1],2) + 
                     math.pow(position1[2] - position2[2],2))                

#finds if a triplet could have an Si atom between them
def dists(positions, dist):

    #if there were not enough close to make a triplet, return none
    if len(positions)<3:
        return[""]
    #if there is a triplet and they are close enough to have a Si, 
    #return the triplet, else return blank
    if len(positions)==3:
        if distance(positions[1],positions[2])<=dist:
            return positions
        else:
            return[""]
    numbers=[]

    if len(positions)==5:
        print(1)
    #if there are more then 2 close enough to have a Si between them,findthe
    # one that could not given the other two
    for i in range(len(positions)):
        numbers.append(0)
    for i in range(1,len(positions)-1):
        for j in range(1,len(positions)-i):
            #if two positions are not close enough, add a counter to both. If they
            # are close enough, remove a counter from both
            if distance(positions[i],positions[i+j])>dist:
                numbers[i]+=1
                numbers[i+j]+=1
            else:
                numbers[i]-=1
                numbers[i+j]-=1

    #removetheonewiththemostcounters
    del positions[numbers.index(max(numbers))]

    #if these still are not close enough to have a triplet between them, return
    # none. If they are close enough, return the new triplet
    if distance(positions[1],positions[2])<=dist:
        return positions
    else:
        return[""]

#finds four membered rings and returns a list of lists of their locaitons
def find_four(opositions, far):
    rings=[[]]
    remov=[]
    #for each oxygen
    for i in range(len(opositions)):
        rings.append([""])
        rings[i]=[opositions[i]]
        #for each oxygen with an x position higher than the current
        for j in range(1,len(opositions)-i):
            #if th exposition is less than the possible distance between two
            # oxygenatoms(variableinclusionradius)
            if abs(opositions[i][0]-opositions[i+j][0])<=far:
                #if the distance between the two oxygens is less than the
                # characteristic distance(variable inclusion radius)
                if distance(opositions[i],opositions[i+j])<=far:
                    rings[i].append(opositions[i+j])
        rem=0
        if len(rings[i])<4:
            rem=1
        elif len(rings[i])>4:
            while len(rings[i])!=4:
                distances=[]
                for k in range(len(rings[i])):
                    tot_len = 0
                    for l in range(1,len(rings[i])-k):
                        tot_len+=distance(rings[i][k],rings[i][k+l])
                    distances.append(tot_len)
                del rings[i][distances.index(max(distances))]
        if len(rings[i])==4:
            distances=[]
            for n in range(len(rings[i])-1):
                for m in range(1,len(rings[i])-n):
                    distances.append(distance(rings[i][n],rings[i][n+m]))
            for n in range(2):
                del distances[distances.index(max(distances))]
            for n in range(4):
                for m in range(1,len(distances)-n):
                    if abs(distances[n]-distances[n+m])>.03:
                        rem=1
        if rem==1:
            remov.insert(0,i)

    for n in range(len(remov)):
        del rings[remov[n]]

    return rings

#finds the area of triangle
def triarea(p1,p2,p3):
    a=distance(p1,p2)
    b=distance(p2,p3)
    c=distance(p1,p3)
    s=(a+b+c)/2
    return math.sqrt(s*(s-a)*(s-b)*(s-c))

def ringarea(corners):
    n=len(corners)
    area=0.0
    for i in range(n):
      j=(i+1)%n
      area += corners[i][0]*corners[j][1]
      area -= corners[j][0]*corners[i][1]
    area=abs(area)/2.0
    return float(area)

#finds if the silicon atom is within a 4 membered ring
def rem4(rings,si):
    deletes=[]
    for i in range(len(rings)):
        triangles=0
        distances=[]
        locations=[]
        for n in range(len(rings[i])-1):
            for m in range(1,len(rings[i])-n):
                distances.append(distance(rings[i][n],rings[i][n+m]))
                locations.append([n,n+m])
        locations.append(len(rings[i]))
        for n in range(2):
            del locations[distances.index(max(distances))]
            del distances[distances.index(max(distances))]
        for n in range(len(locations)):
            triangles += triarea(rings[i][locations[n][0]],
                                 rings[i][locations[n][1]],si)
        if ringarea(rings[i])==triangles:
            return"n"
    return"y"

#finds the position of a Si given a triplet of oxygen
def si_finder(opositions):

    #characteristic distance
    dist=1.6*math.pow(10,-1)

    #sets up the translation to happen around a basepoint(the first point in
    # the positions)
    trans=[[0,0,0],[opositions[1][0]-opositions[0][0], 
            opositions[1][1]-opositions[0][1],
            opositions[1][2]-opositions[0][2]],
    [opositions[2][0]-opositions[0][0],opositions[2][1]-opositions[0][1],
     opositions[2][2]-opositions[0][2]]]

    #finds vector perpendicular to the plane of the three points
    v=numpy.matrix([numpy.linalg.det([[trans[1][1],trans[2][1]],
                                      [trans[1][2],trans[2][2]]]),
    numpy.linalg.det([[trans[1][0],trans[2][0]],[trans[1][2],
                       trans[2][2]]]),
    numpy.linalg.det([[trans[1][0],trans[2][0]],[trans[1][1],trans[2][1]]])])

    #sets up first rotation matrix about the x axis
    theta=math.atan2(v.item(1),v.item(2))
    xmatr=numpy.matrix([[1,0,0],[0,math.cos(theta),-math.sin(theta)],
                        [0,math.sin(theta),math.cos(theta)]])
    trans1=numpy.matrix(trans)
    rot1=numpy.matrix.dot(trans1,xmatr)
    v1=numpy.matrix.dot(v,xmatr)

    #second rotation matrix about the y axis
    rho=math.atan2(v1.item(0),v1.item(2))
    ymatr=numpy.matrix([[math.cos(rho),0,math.sin(rho)],[0,1,0],
                         [-math.sin(rho),0,math.cos(rho)]])
    rot2=numpy.matrix.dot(rot1,ymatr)

    #should be in the xy plane now. Have to rotate such that two points 
    # are on the x axis
    alph=math.atan2(rot2.item(4),rot2.item(3))
    bet=math.atan2(rot2.item(7),rot2.item(6))
    r1=math.sqrt(math.pow(rot2.item(3),2)+math.pow(rot2.item(4),2))
    r2=math.sqrt(math.pow(rot2.item(6),2)+math.pow(rot2.item(7),2))
    rot3=numpy.matrix([[rot2.item(0),rot2.item(1),rot2.item(2)],[r1,0,0],
                        [r2*math.cos(bet-alph),r2*math.sin(bet-alph),0]])
    x=r1/2
    y=r2*(1-math.cos(bet-alph))/(2*math.sin(bet-alph))
    z=math.sqrt(abs(math.pow(dist,2)-math.pow(x,2)-math.pow(y,2)))
    si_pos=numpy.matrix([x,y,z])

    #rotate back to originial position
    init=math.atan2(si_pos.item(1),si_pos.item(0))
    r=math.sqrt(math.pow(si_pos.item(0),2)+math.pow(si_pos.item(1),2))
    x=r*math.cos(init+alph)
    y=r*math.sin(init+alph)
    si_pos=numpy.matrix([x,y,z])

    #undo second rotation matrix
    iymatr=numpy.linalg.inv(ymatr)
    si_pos=numpy.matrix.dot(si_pos,iymatr)

    #undo first rotation matrix
    ixmatr=numpy.linalg.inv(xmatr)
    si_pos=numpy.matrix.dot(si_pos,ixmatr)

    #translate back so there is no point at the origin
    si_pos=[si_pos.item(0)+opositions[0][0],si_pos.item(1)+opositions[0][1], 
            si_pos.item(2)+opositions[0][2]]

    return si_pos
    
#locates all possiblee triplets
def o_locator(opositions):
    
    dist = 1.6*math.pow(10,-1)
    
    #assumed oxygens are ordered by increasing x values
    #used to collect all the found oxygens close enough to have a single Si
    # between them
    found=[[""]]
    #for each oxygen
    for i in range(len(opositions)):
        found[i]=[opositions[i]]
        #for each oxygen with an x position higher than the current
        for j in range(1,len(opositions)-i):
            #if the x position is less than the possible distance between two
            # oxygenatoms(variableinclusionradius)
            if abs(opositions[i][0]-opositions[i+j][0])<=3.45*math.pow(10,-1):
            #if the distance between the two oxygens is less than the
                # characteristic distance(variable inclusion radius)
                if distance(opositions[i],opositions[i+j])<=3.45*math.pow(10,-1):
                    found[i].append(opositions[i+j])
        found.append([""])

    #removes last appended empty list
    del found[len(found)-1]

    #remove all those too far apart using dist function (variable inclusion
    # radius)
    for n in range(len(found)):
        found[n]=dists(found[n], .345)


    #createanarrayforpositionstoremove
    remov=[]
    #for all atoms with found oxygens
    for n in range(len(found)):
        #add empties to a list for removal
        if found[n]==[""]:
            remov.insert(0,n)

    #remove those in the remove list
    for m in range(len(remov)):
        del found[remov[m]]

    #return the list of those oxygen that have a possible Si between them
    return found


def main():

    #convert data in file into floats and append to a position list
    with open("OfC Positions 120106_008 Python Output_2.txt") as f:
        content=f.readline()

    string = ""

    locations = []
    
    for i in range(len(content)):
        if content[i] == " ":
            locations.append(float(string))
            string = ""
        else:
            string += content[i]

    locations.append(float(string))

    positions=[[""]]

    for i in range(len(locations)):
        if i%3==0:
            positions[i/3]=[locations[i]]
            positions.append("")
        else:
            positions[i/3].append(locations[i])

    del positions[len(positions)-1]
     
    #find triplets
    triples=o_locator(positions)
 
     #find Si positions
    si_locations=[]
    for j in range(len(triples)):
         si_locations.append(si_finder(triples[j]))
 
     #rings = find four(positions,.35)
 
    delete=[]
 
     #for i in range(len(silocations)):
         #if rem4(rings,silocations[i])=="n":
             #delete.append(0,i)
 
    for i in range(len(delete)):
        del si_locations[delete[i]]
    
    print(si_locations)
    
    #write Si positions to an outfile
    out = open("test.txt", "w")
    out.write(str(si_locations))
    out.write("\n")
 
if __name__=="__main__": main()
