import numpy as np

def find2d(coords,coord):
        #Find the nearest value from a 2D coordinates.
        return np.argmin(abs(np.sum(coords-coord,axis=1)))
    
def find(coords,coord):
    #Find nearest value `coord` from vector `coords`
    return np.argmin(abs(coords-coord))

def slope(X1,X2,Y1,Y2):
    #Compute slope
    return (Y2-Y1)/(X2-X1)

def slope2angle(slopes):
    return np.arctan(slopes)