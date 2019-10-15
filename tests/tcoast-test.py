import numpy as np
from tcoasts import *
import pytest
import os
import gsw

package_path=tcoasts.__file__ .split('tcoasts.py')[0]

folder=os.path.join(package_path,'../examples/data/')
contourfile='GoM_coastLine.xy'

@pytest.mark.ttcoasts
def test_perpendicularity():
    dotp=perpendicularity(method='smooth')
    print(dotp)
    assert (dotp<1e-10).all() 

def perpendicularity(method):
    '''
    Test suit 1: Ensure that perpendicular vectors are perpendicular.
    '''
    distance=np.arange(0,200,50)

    tac=tcoasts.TransportAlongCoast(folder,[-89.75,21.3],contourfile,distance)
    
    vec_dict=tac.perp2coast(method=method)
    
    perp_vect=np.vstack((tac.x_perp,tac.y_perp)).T
    norm_vect=np.vstack((tac.x_norm,tac.y_norm)).T

    dotp=np.array([np.dot(perp_vect[ii],norm_vect[ii]) for ii in range(len(perp_vect))])
    return dotp


@pytest.mark.parametrize(('method'), [
    ('smooth'),('byseg'),('local'),('ext')
])

@pytest.mark.ttcoasts
def test_perpendicularity_methods(method):
    dotp=perpendicularity(method=method)
    assert (dotp<1e-10).all() 

########## Add test to know perpendicular vel vector and the normal projection

########## Add test if file exists 
    
@pytest.mark.parametrize(('n'), [
    (10),(1)
])

@pytest.mark.ttcoasts
def test_distance(n):
    dist,length=testdist(n)
    assert abs(1-((length)/dist)).max() < 1e-3

def testdist(n):
    distance=np.arange(0,200,n)
    tac=tcoasts.TransportAlongCoast(folder,[-89.75,21.3],contourfile,distance)
    locations = tac.perploc()
    x_norm,y_norm,x_perp,y_perp,x_perp_all,y_perp_all=tac.vertor_perp()
    dist=np.zeros(len(locations))
    for ii in range(len(dist)):
        dist[ii]=gsw.distance([x_perp_all[ii][0],x_perp_all[ii][-1]],[y_perp_all[ii][0],y_perp_all[ii][-1]])
    return dist,tac.length*1e3

@pytest.mark.parametrize(('transect'), [
    (None),(0)
])

@pytest.mark.ttcoasts
def test_plot(transect):
    testplot(transect)

def testplot(transect):
    distance=np.arange(0,200,2)
    tac=tcoasts.TransportAlongCoast(folder,[-89.75,21.3],contourfile,distance)
    tac.plotperp_vect(transect=transect)

@pytest.mark.ttcoasts
def test_transport():
    testtransport()

def testtransport():
    distance=np.arange(0,200,10)
    tac=tcoasts.TransportAlongCoast(folder,[-89.75,21.3],contourfile,distance)
    tac.inter2vector(ufiles='U.nc',vfiles='V.nc',save=False)
    tac.compute_transport()

@pytest.mark.ttcoasts
def test_interp_decorator():
    for ii in range(0,2):
        testinterpdecorator()
    fileexists=os.path.isfile('./tmp_interp_transects.nc')
    os.remove('./tmp_interp_transects.nc')
    assert fileexists
    

def testinterpdecorator():
    distance=np.arange(0,200,10)
    tac=tcoasts.TransportAlongCoast(folder,[-89.75,21.3],contourfile,distance)
    tac.inter2vector(ufiles='U.nc',vfiles='V.nc',save=True)