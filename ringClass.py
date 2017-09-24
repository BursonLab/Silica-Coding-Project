## Create the general ring class

# IMPORTS
import graphics

#win = GraphWin('Test Remove', 500, 500)

class Ring:

    """A button is a labeled curcle in a window.
    It is activated or deactivated with the activate()
    and deactivate() methods. The clicked(p) method
    returns true if the button is active and p is inside it."""

    def __init__(self, win, center, radius, color, label):
        """ Creates a circular ring, eg:
        qb = Button(myWin, Point(30,25), 20, 'blue', 'Quit') """ 

        self.color = color
        self.outline = 'black'
        
        self.center = center
        self.x,self.y = self.center.getX(), self.center.getY()
        self.radius = radius

        
        self.circle = graphics.Circle(self.center, radius)
        self.circle.draw(win)
        
        self.label = label

        self.text = graphics.Text(self.center, self.label)
        self.text.setSize(8)
        self.text.draw(win)

        self.activate()
        
        
    def inside(self, p):
        """RETURNS ture if click in inside button"""
        return self.active and ((self.x - p.getX())**2 + (self.y-p.getY())**2)**(1/2) < self.radius

    def clicked(self, p):
        """ RETURNS true if button active and p(clickpoint -> win.getMouse()) is inside"""
        if self.inside(p):
            self.circle.setOutline('yellow')
            self.outline = 'yellow'
            return True
        else:
            return False

    def getLabel(self):
        """RETURNS the label string of this button."""
        return self.label.getText()

    def activate(self):
        """Sets this button to 'active'."""
        self.circle.setFill(self.color)
        self.circle.setOutline('black')
        self.outline = 'black'
        self.circle.setWidth(2)
        self.active = 1
        
    def deactivate(self):
        """Sets this button to 'inactive'."""
        self.circle.setFill('lightgrey')
        self.circle.setOutline('black')
        self.outline = 'black'
        self.circle.setWidth(1)
        self.active = 0
        
    def getOutlineColor(self):
        return self.outline
    
    def removeRing(self):
        self.circle.undraw()
        self.text.undraw()
        
    def getCenter(self):
        return self.center
    
#    def getLabel(self):
#        return self.label
        
        
##def main(win):
##    ring = Ring(win, Point(250, 250), 50, 'green', '4MR')
##    if ring.clicked(win.getMouse()):
##        ring.remove()
##                    
##main(win)

    

