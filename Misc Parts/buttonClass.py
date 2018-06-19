## Create the general button class

# IMPORTS
import graphics

class Button:

    """A button is a labeled rectangle in a window.
    It is activated or deactivated with the activate()
    and deactivate() methods. The clicked(p) method
    returns true if the button is active and p is inside it."""

    def __init__(self, win, center, width, height, color, label):
        """ Creates a rectangular button, eg:
        qb = Button(myWin, Point(30,25), 20, 10, 'Quit') """ 

        self.color = color
        self.outline = 'black'
        
        w,h = width/2.0, height/2.0
        x,y = center.getX(), center.getY()
        self.xmax, self.xmin = x+w, x-w
        self.ymax, self.ymin = y+h, y-h
        p1 = graphics.Point(self.xmin, self.ymin)
        p2 = graphics.Point(self.xmax, self.ymax)
        
        self.rect = graphics.Rectangle(p1,p2)
        self.rect.draw(win)
        
        self.label = graphics.Text(center, label)
        
        self.label.draw(win)
        self.activate()
        
        
    def inside(self, p):
        """RETURNS ture if click in inside button"""
        return self.active and \
               self.xmin <= p.getX() <= self.xmax and \
               self.ymin <= p.getY() <= self.ymax
        

    def clicked(self, p):
        """ RETURNS true if button active and p(clickpoint -> win.getMouse()) is inside"""
        if self.inside(p):
            self.rect.setOutline('yellow')
            self.outline = 'yellow'
            return True
        else:
            return False

    def getLabel(self):
        """RETURNS the label string of this button."""
        return self.label.getText()

    def activate(self):
        """Sets this button to 'active'."""
        self.rect.setFill(self.color)
        self.rect.setOutline('black')
        self.outline = 'black'
        self.rect.setWidth(2)
        self.active = 1
        
    def deactivate(self):
        """Sets this button to 'inactive'."""
        self.rect.setFill('lightgrey')
        self.rect.setOutline('black')
        self.outline = 'black'
        self.rect.setWidth(1)
        self.active = 0
        
    def getOutlineColor(self):
        """Returns outline color"""
        return self.outline
        
    def setBoldText(self):
        """Set the button label bold"""
        self.label.setStyle('bold')
        
    def setSizeText(self, size):
        """Set the size of the label"""
        self.label.setSize(size)
        

    
