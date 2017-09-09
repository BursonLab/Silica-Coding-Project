"""***************************************************************************************\
PROGRAM NAME:         ringAnalyzer.py

PROGRAM DESCRPTION:   Program to place rings/atoms on an image and returns the center 
                      position of each placed object in pixels and/or nm
******************************************************************************************"""

#IMPORTS
from graphics import *
from Image import PIL
from buttonClass import Button
from ringClass import Ring
import matplotlib.pylab as plt
import numpy as np


#Set up main windown
win = GraphWin('Ring Analyzer', 1500, 1000)

def createBoard(win):
    """Creates program boad & instructions"""
    
    #Create line diving instructions, toolbar, and image
    divideLine1 = Line(Point(498, 0), Point(498, 1000))
    divideLine1.setWidth(2)
    divideLine1.draw(win)
    
    divideLine2 = Line(Point(0, 600), Point(498, 600))
    divideLine2.setWidth(2)
    divideLine2.draw(win)
    
    instructionLable = Text(Point(250, 25), 'Instructions')
    instructionLable.setSize(20)
    instructionLable.setStyle('italic')
    instructionLable.draw(win)
    
    toolLabel = Text(Point(250, 625), 'Toolbar')
    toolLabel.setSize(20)
    toolLabel.setStyle('italic')
    toolLabel.draw(win)
    
    imageLabel = Text(Point(1000, 500), 'Upload Image Here')
    imageLabel.setSize(20)
    imageLabel.setStyle('italic')
    imageLabel.draw(win)
    
    #INSTRUCTIONS
    instruction1 = Text(Point(250, 75), "1.) Enter the name of the image file, then click ENTER. \
    \n Note that the file must be saved as .png \n (ex. filename.png)")
    instruction1.draw(win)

    entry1 = Entry(Point(250, 125), 25)
    entry1.draw(win)
    
    enterButton1 = Button(win, Point(425, 125), 80, 20, '#52E643', 'ENTER')
    
    #Click enterButton, open the image, deactivate image
    loop = True
    while loop:
        if enterButton1.clicked(win.getMouse()):
            imageName = entry1.getText()
            #Check user imput for name of file
            if imageName != '' and imageName[-4:] == '.png':
                loop = False
            else:    
                errorMessage('Please enter a PNG image file. \n (ex. filename.png)')
                enterButton1.activate()
                
    #Upload image into window
    myImage = Image(Point(1000, 500), imageName)
    myImage.draw(win)
    enterButton1.deactivate()
                    
    instruction2 = Text(Point(250, 190), "2.) Click on the buttons below to place an atom or ring center. \
    \n Double click to place the object on the image. \
    \n Then click FINISH when all of the objects have been placed.")
    instruction2.draw(win)
    
    #Create TOOLBAR
    siButton = Button(win, Point(300, 700), 40, 25, '#1DAA43', 'Si')
    oButton = Button(win, Point(400, 700), 40, 25, '#FF1D0F', 'O')
    
    generalRing = Button(win, Point(100, 680), 140, 34, '#F8A61E', 'General Ring')
    button4MR = Button(win, Point(100, 750), 140, 22, '#9E4BF6', '4 Membered Ring')
    button5MR = Button(win, Point(100, 785), 140, 22, '#4BB0F6', '5 Membered Ring')
    button6MR = Button(win, Point(100, 820), 140, 22, '#43C47F', '6 Membered Ring')
    button7MR = Button(win, Point(100, 855), 140, 22, '#F9DE3E', '7 Membered Ring')
    button8MR = Button(win, Point(100, 890), 140, 22, '#F94C3E', '8 Membered Ring')
    button9MR = Button(win, Point(100, 925), 140, 22, '#F726E8', '9 Membered Ring')
   
    
    removeButton = Button(win, Point(350, 845), 100, 25, '#F96A61', 'REMOVE')
    doneButton = Button(win, Point(350, 800), 100, 25, '#41EFC9', 'DONE')
    finishButton = Button(win, Point(350, 900), 120, 30, '#D4FF33', 'FINISH')
    
    
    return [siButton, oButton, generalRing, button4MR, button5MR, button6MR, button7MR, button8MR, \
           button9MR, removeButton, finishButton, doneButton]


def checkDone(mousePress, doneButton):
    """Check if DONE button is clicked. If yes, return TRUE"""
    
    if doneButton.clicked(mousePress):
        return True
    else:
        return False
    
    
def resetButtons(buttonLst):
    """Activates buttons in buttonLst"""
    
    for i in buttonLst:
        i.activate()
        
        
def placeCenter(win, clickButton, doneButton, centerPtLst, genRingLst, buttonsLst, color, label):
    """Places a center when the corresponding button is pressed"""

    #Deactivate other buttons (not DONE or the selected button)
    deactButtons = buttonsLst[:-1]
    deactButtons.remove(clickButton)
    for i in deactButtons:
        i.deactivate()
                
    while clickButton.getOutlineColor() == 'yellow':
        mousePress = win.getMouse()
        if checkDone(mousePress, doneButton) == False:
            
            centerPtLst.append(win.getMouse())
        
            ring = Ring(win, centerPtLst[-1], 10, color, label)
            genRingLst.append(ring)
            
        elif checkDone(mousePress, doneButton) == True:
            resetButtons(buttonsLst)
                    
        
def convertCenter2nm(pixelCenter, conversionFactor):
    """Converts center point in pixels to nanometers given facor in nm/pixels"""
    
    x, y = pixelCenter.getX(), pixelCenter.getY()
    return (x * conversionFactor, y * conversionFactor)


def plotDistribution(labelLst, colorLst, centerPointsLst):
    """Plots the distribution of atoms and ring sizes
    Note:  centerPointsLst is a nested list of all center points"""

    datas = []
    for i in range(len(labelLst)):
        datas.append({'label': labelLst[i], 'color': colorLst[i], 'height': len(centerPointsLst[i])})  
    
    i = 0
    for data in datas:
        plt.bar(i, data['height'],align='center',color=data['color'])
        i += 1
    
    labels = [data['label'] for data in datas]
    pos = [i for i in range(len(datas)) ]
    plt.xticks(pos, labels)
    plt.xlabel('Atom / Ring Size')
    plt.ylabel('Number of Atoms / Rings')
    plt.title('Distribution of Atoms / Ring Sizes')
    plt.show()
    
    
def errorMessage(message):
    """Displays an error message in a new window"""
    
    win = GraphWin('Error Message', 350, 175)
    win.setBackground('#F75454')
    
    
    error = Text(Point(win.getWidth()/2, win.getHeight()/2 - 50), 'Error')
    error.setSize(30)
    error.setStyle('bold')
    error.draw(win) 
    
    text = Text(Point(win.getWidth()/2, win.getHeight()/2), message)
    text.setFace('times roman')
    text.draw(win)
    
    #Press OK button to close the window
    okButton = Button(win, Point(win.getWidth()/2, win.getHeight()/2 + 55), 100, 30, 'lightgrey', 'OK')
    loop = True
    while loop:
        mousePress = win.getMouse()
        if okButton.clicked(mousePress):
            win.close()
            loop = False
            
            
def alertMessage(message):
    """Displays an error message in a new window"""
    win = GraphWin('Alert Message', 300, 100)
    win.setBackground('#6FF3CC')
    
    text = Text(Point(win.getWidth()/2, win.getHeight()/2 - 20), message)
    text.setSize(11)
    text.setFace('times roman')
    text.draw(win)
    
    
    okButton = Button(win, Point(win.getWidth()/2 - 70, win.getHeight()/2 + 30), 80, 25, 'lightgrey', 'OK')
    okButton.setSizeText(10)
    cancelButton = Button(win, Point(win.getWidth()/2 + 70, win.getHeight()/2 + 30), 80, 25, 'lightgrey', 'CANCEL')
    cancelButton.setSizeText(10)
    
    #Press OK button to continue with function or cancel to quit alert window
    loop = True
    while loop:
        mousePress = win.getMouse()
        if okButton.clicked(mousePress):
            win.close()
            loop = False
            return True
        
        elif cancelButton.clicked(mousePress):
            win.close()
            loop = False
            return False
            
            
def checkDecimal(string):
    """Checks if the string has a decimal place"""
    
    for i in string:
        if i == '.':
            return True
    return False
    
    
def writeXYZ(filename, atomLst, centerLst):
    """Exports center poitision to an xyz file"""
    
    #Convert nm -> angstroms    
    for i in range(len(centerLst)):
        for j in range(len(centerLst[i])):
            x, y = centerLst[i][j]
            x = x*10
            y = y *10
            centerLst[i][j] = x, y

   
    n = 0   #Varibale to count number of atoms
    #Convert points from touple -> list
    #Set z coordinate to 0 
    for i in range(len(centerLst)):
        for j in range(len(centerLst[i])):
            centerLst[i][j] = list(centerLst[i][j])
            centerLst[i][j].append(0.00)
            n += 1  #Update the count
            
            #Round each point to 5 decimal places
            for k in range(len(centerLst[i][j])):
                centerLst[i][j][k] = round(centerLst[i][j][k], 5)
                
                #Make each point the same length
                while len(str(centerLst[i][j][k])) != 9:
                    centerLst[i][j][k] = str(centerLst[i][j][k]) + '0'

    #Format the xyz file
    with open(filename, 'w') as file:
   
        file.write(str(n) + '\n')
        for i in range(len(centerLst)):
            
            for j in range(len(centerLst[i])):
                file.write('\n' + str(atomLst[i]))
                file.write('    ' + str(centerLst[i][j][0]))
                file.write('    ' + str(centerLst[i][j][1]))
                file.write('    ' + str(centerLst[i][j][2]))





def main(win):
    """The main code of the program"""
    
    #Create the board & buttons
    buttons = createBoard(win)
    
    # Make a nested list for ring centers    
    centerLst = []
    for i in range(9):
        centerLst.append([])    
    #centerLst = [[Si Centers], [O], [General Ring], [4 Membered], [5 Membered], [6 Membered],
    #              [7 Membered], [8 Membered], [9 Membered]]
        
        
    # Define Lists to be used later
    colorLst = ['#1DAA43', '#FF1D0F', '#F8A61E', '#9E4BF6', '#4BB0F6', \
                '#43C47F', '#F9DE3E', '#F94C3E', '#F726E8']            #Colors for rings
    labelLst = ['Si', 'O', 'GR', '4', '5', '6', '7', '8', '9']      #Labels for rings
    ringLst = []                                                 #List to store ring objects
    
    
    #Check if a button is clicked & place rings until user clicks FINISH
    loop1 = True
    while loop1:
        mousePress = win.getMouse()
        
        #Place rings when button is clicked
        for i in range(9):
            if buttons[i].clicked(mousePress):
                placeCenter(win, buttons[i], buttons[11], centerLst[i], \
                            ringLst, buttons, colorLst[i], labelLst[i])
                
        #REMOVE atoms/rings when button is clicked
        if buttons[9].clicked(mousePress):
            removeRingsLst = []
            
            #Deactivate all buttons except REMOVE and DONE
            for i in range(9):
                buttons[i].deactivate()
            buttons[10].deactivate()
                
            #Select rings to remove
            while buttons[9].getOutlineColor() == 'yellow':
                mousePress = win.getMouse()
                for i in ringLst:
                    if i.clicked(mousePress):
                        i.removeRing()
                        ringLst.remove(i)
                        
                        #Remove center from corresponding list
                        center = i.getCenter()
                        for j in range(len(labelLst)):
                            if i.getLabel() == labelLst[j]:
                                centerLst[j].remove(center)
                            
                if checkDone(mousePress, buttons[11]) == True:
                    resetButtons(buttons)
                    
        #Input dimensions of photo when FINISH is clicked
        elif buttons[10].clicked(mousePress):
            #Exist fist while loop
            loop1 = False
                
            #Deactivate all buttons
            for i in buttons:
                i.deactivate()
            
            #Find pixel -> nm conversion factor           
            #Enter Pixel dimensions
            instruction3 = Text(Point(250, 265), "3.) Enter the size of the image file in pixels. \
            \n Note that the entry must include a decimal point.\n(ex. 900.0 x 900.0 pixels)")
            instruction3.draw(win)
             
            #Enter nm dimensions
            instruction4 = Text(Point(250, 360), "4.) Enter the size of the image file in nanometers (nm). \
            \n Note that the entry must include a decimal point.\n(ex. 5.0 x 5.0 nm)")
            instruction4.draw(win)
                            
            #Set up entry boxes for user
            nmBoxH = Entry(Point(282, 405), 7)
            nmBoxH.draw(win)
            
            text2 = Text(Point(237.5, 405), "X")
            text2.draw(win)
            
            txtLabel2 = Text(Point(330, 405), "nm")
            txtLabel2.draw(win)
        
            nmBoxW = Entry(Point(194, 405), 7)
            nmBoxW.draw(win)
                
            pixBoxH = Entry(Point(282, 310), 7)
            pixBoxH.draw(win)
        
            text = Text(Point(237.5, 310), "X")
            text.draw(win)
        
            txtLabel = Text(Point(340, 310), "pixels")
            txtLabel.draw(win)
        
            pixBoxW = Entry(Point(194, 310), 7)
            pixBoxW.draw(win)

            enterButton2 = Button(win, Point(425, 405), 85, 20, '#52E643', 'ENTER')
        
            #Program dimensions enter button
            loop2 = True
            while loop2:
                mousePress = win.getMouse()

                if enterButton2.clicked(mousePress):
                    #enterButton2.deactivate()
                    
                    #Sore values of dimensions
                    nmHeight = nmBoxH.getText()
                    nmWidth = nmBoxW.getText()
                    pixHeight = pixBoxH.getText()
                    pixWidth = pixBoxW.getText()
                    
                    #Check the user input values are decimals
                    if checkDecimal(nmHeight) and checkDecimal(nmWidth) and checkDecimal(pixHeight) \
                    and checkDecimal(pixWidth):
                
                        #Calculate average conversion factor
                        convtFact = (((float(nmWidth) / float(pixWidth)) + \
                        (float(nmHeight) / float(pixHeight))) / 2)

                        #Finish the program after user has input dimensions
                        instruction5 = Text(Point(250, 450), "5.) Click on the buttons below to finish the program.")
                        instruction5.draw(win)
                        
                        enterButton2.deactivate()
                    
                        #Exit second while loop
                        loop2 = False
                    
                    else:
                        errorMessage('Please enter dimensions as decimals. \n (ex. 900.0 x 900.0 pixels)')
                        enterButton2.activate()
                    

    #Create final buttons
    exportButton = Button(win, Point(100, 520), 100, 40, '#FFEE33', 'Export to \n xyz File')
    exportButton.setBoldText()
    distButton = Button(win, Point(250, 520), 150, 40, '#33D6FF', 'Show Distribution \n Graph')
    distButton.setBoldText()
    quitButton = Button(win, Point(400, 520), 100, 40, '#FF3369', 'QUIT')
    quitButton.setBoldText()
    
    
    #List of labels
    textLst = ['Si', 'O ', 'GR', '4M', '5M', '6M', '7M', '8M', '9M']
    
    #Print lists in pixels
    #for i in range(len(textLst)):
    #    print(textLst[i] + ' Centers (pixels) = ' + str(centerLst[i]))
        
    #Convert list of centers in pixels -> nm
    for i in range(len(centerLst)):
        for j in range(len(centerLst[i])):
            centerLst[i][j] = convertCenter2nm(centerLst[i][j], convtFact)
        
    #Print lists in nm
    #for i in range(len(textLst)):
    #    print(textLst[i] + 'Centers (nm) = ' + str(centerLst[i]))
          
    #Program last three buttons
    loop3 = True
    while loop3:
        mousePress = win.getMouse()
        if exportButton.clicked(mousePress):
            #Export to xyz file (cords in Angstroms)
            if alertMessage('Export coordinates to xyz file \n saved as ringAnalyzerCoordinates.txt'):
                writeXYZ('ringAnalyzerCoordinates.txt', textLst, centerLst)
                exportButton.deactivate()
            else:
                exportButton.activate()
                            
        elif distButton.clicked(mousePress):
            #PLOT atom / ring distribution
            plotDistribution(textLst, colorLst, centerLst)
            distButton.deactivate()
                            
                            
        elif quitButton.clicked(mousePress):
            if alertMessage('Click OK to exit the program \n (The image file will not be saved)'):
                win.close()
                loop3 = False
            else:
                quitButton.activate()
            

main(win)