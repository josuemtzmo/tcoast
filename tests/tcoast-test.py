import numpy as np
import tcoasts

def perpendicularity():
    ### Function needs to be updated
    '''
    Test suit 1: Ensure that perpendicular vectors are perpendicular.
    '''
    folder='/g/data3/hh5/tmp/jm_eddy_lagrangian/mitgcm/run/'
    contourfile='./input_data/GoMCoastLine_nolagoon.xy'
    
    tac=TransportAlongCoast(folder,[-94,18],contourfile)
    
    x_norm,y_norm,x_perp,y_perp=tac.perp2coast(method='smooth')
    
    a=np.array([x_norm[ii]-tac.coastline[locations[ii],0],y_norm[ii]-tac.coastline[locations[ii],1]]).T
    b=np.array([x_perp[ii]-tac.coastline[locations[ii],0],y_perp[ii]-tac.coastline[locations[ii],1]]).T
    # Make sure this number is always smaller than 1e-10 (i.e. zero)
    np.dot(a,b)
    
### add test which will test different methods.
    
def testdist():
    ### Function needs to be updated
    folder='/g/data3/hh5/tmp/jm_eddy_lagrangian/mitgcm/run/'
    contourfile='./input_data/GoMCoastLine_nolagoon.xy'
    tac=TransportAlongCoast(folder,[-94,18],contourfile)
    locations=tac.perploc()
    x_norm,y_norm,x_perp,y_perp=tac.perp2coast(method='smooth')
    dist=np.zeros(len(locations))
    for ii in range(len(locations)):
        dist[ii]=gsw.distance([x_perp[ii][0],x_perp[ii][-1]],[y_perp[ii][0],y_perp[ii][-1]])

    if (dist > 1e3*(tac.length+tac.length*1e-4)).any():
        raise ValueError('Error, the perpendicular distance has not been computed correctly.')
       