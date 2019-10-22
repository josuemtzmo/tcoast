import gsw
import xarray as xr
import subprocess
import numpy as np
import os
import pylab as plt
# Import utils and decorators
from tcoasts.utils.utils import *
from tcoasts.utils.decorators import _file_exists


class TransportAlongCoast(object):
    '''
    
    '''
    def __init__(self,path,initpos,contour_file,distance=np.arange(-400,400,100),length=100):
        self.path     = path # Data path
        self.initpos  = initpos # Init location lat lon coordinates.
        self.dist     = distance # Units kilometers
        self.length   = length # Units kilometers
        self.n        = 4 # self.length/self.n corresponds to the segments 
                          # on the perpendicular vector.
        self.contour_file = contour_file # Contour filename.
        self.tmpfile= 'tmp_interp_transects.nc' # Temporal file to store interpolated fields.
        # Load required data
        self.extract_contour()
        self.initindex=find2d(self.coastline,initpos)
    
    def extract_contour(self):
        # Load contour file.
        if './' not in self.contour_file:
            self.contour_file=os.path.join(self.path,self.contour_file)
        if os.path.isfile(self.contour_file):
            self.coastline=np.loadtxt(self.contour_file)
        else:
            raise ValueError('''
                                Make sure the file path is correct. 
                                The path should be relative to the location of 
                                the running script, or relative to self.path.
                            ''')
    
    def coords2dist(self,lonlat,p=0):
        '''
        This function follows the GSW computation.
        '''
        distance=gsw.distance(lonlat[:,0],lonlat[:,1],p)
        return distance
    
    def distancefrominit(self):
        '''
        The distance definition points positive to the east.
        '''
        if self.initindex != 0:
            # Compute cumulative distance to right of index [location,location]
            postinit=np.cumsum(self.coords2dist(self.coastline[self.initindex:]))
            # Compute cumulative distance to left of index [location,location]
            neginit=-1*np.cumsum(np.flipud(self.coords2dist(self.coastline[:self.initindex])))
            # Join cumulative distances.
            cumdistance=np.hstack((np.flipud(neginit),postinit))
        else:
            # Compute cumulative distance starting from the index [0,0]
            cumdistance=np.cumsum(self.coords2dist(self.coastline))
        return cumdistance

    def perploc(self):
        #Find the user defined locations for the perpendicular vectors.
        dist_coast=self.distancefrominit()
        index_perp=[find(dist_coast,dis*1000) for dis in self.dist]
        return index_perp
    
    def perp2coast(self,method='smooth',x=10):
        '''
        Input:
        
            method: [ mean ]
                smooth  - computes the mean over X number of slopes and 
                          projects the perpendicular vector
                byseg   - computes the mean over each segment of the slope
                local   - computes the the perpendicular vector using the 2 
                          adjacent locations
                ext     - computes the perpendicular vector using the slope at 
                          x cells to the left and right of the desired 
                          perpendicular location.
        
        '''
        index_perp=self.perploc()
        # Method to find the location perpendicular vector.
            
        if method =='local' and method =='ext':
            # Compute slope from adjacent locations [loc-x,loc+x]
            if method=='local':
                x=1
            slopes=np.array([slope(self.coastline[ii-x,0],self.coastline[ii+x,0],
                                        self.coastline[ii-x,1],self.coastline[ii+x,1]) 
                             for ii in index_perp])
            
        elif method == 'smooth':
            # Compute average slope from all the indexes contained between locations [loc-x,loc+x] 
            slopes=np.array([np.mean([slope(self.coastline[ii-xx,0],self.coastline[ii+xx,0],
                                        self.coastline[ii-xx,1],self.coastline[ii+xx,1]) 
                                      for xx in range(1,x)]) 
                             for ii in index_perp])
        
        else:
            # Compute average slope from all segments from [loc-x,loc-x+(2x-1)]
            slopes=np.array([np.mean([slope(self.coastline[ii-x,0],self.coastline[ii-x+xx,0],
                                        self.coastline[ii-x,1],self.coastline[ii-x+xx,1]) 
                                      for xx in range(1,(2*x-1))]) 
                             for ii in index_perp])
        
        #Compute angles from slopes
        angles=slope2angle(slopes)
        #Shift angles to be perpendicular
        perp_angle=angles+(np.pi/2)
        #Normal vector
        self.x_norm = np.squeeze(np.cos(angles))
        self.y_norm = np.squeeze(np.sin(angles))
        
        #Perpendicualar vector
        self.x_perp = np.squeeze(np.cos(perp_angle))
        self.y_perp = np.squeeze(np.sin(perp_angle))
        
        # Return dictionary containing vector information
        return {'Nvector':{'x':self.x_norm,'y':self.x_norm,'angle':angles,'slope':slopes},
                'Pvector':{'x':self.x_perp,'y':self.y_perp,'angles':perp_angle,'slope':-1/slopes}}
    
    def perpvecdist(self,index_perp,perp_angle):
        #compute distances to scale perpendicular vectors.
        ### Note this will produce an error of 1e-4.
        x=np.array([[self.coastline[index_perp][ii,0],
                     np.cos(perp_angle[ii])+self.coastline[index_perp][ii,0]] 
                     for ii in range(len(index_perp))])
        y=np.array([[self.coastline[index_perp][ii,1],
                     np.sin(perp_angle[ii])+self.coastline[index_perp][ii,1]] 
                     for ii in range(len(index_perp))])
        distances = gsw.distance(x,y)
        return distances
    
    # _file_exists will test if the tmporal file containing the interpolated 
    # data exits. If file exists it will load the contents, otherwise, it will 
    # interpolate the data.
    @_file_exists
    def inter2vector(self,ufiles='U.*.nc',vfiles='V.*.nc',tracerfile=None,dataset=None,save=True,**kwargs):
        '''
        **kwargs inter2vector supports the all the  kwargs of xr.open_mfdataset.
        '''
        # xr load parameters
        xr_openmf_defaults={}
        if '*' in ufiles and '*' in vfiles:
            xr_openmf_defaults = {'concat_dim':'time','parallel':True,'combine':'nested'}
            xr_openmf_defaults.update(kwargs)

        print('Opening velocity files')
        if dataset != None:
            # Load data.
            u = self.loaddata(file=ufiles,var='U',dataset=dataset,**xr_openmf_defaults)
            v = self.loaddata(file=vfiles,var='V',dataset=dataset,**xr_openmf_defaults)
        else:
            u = dataset.U
            v = dataset.V
        # Make sure the shape of the velocity fields are the same.
        if u.shape != v.shape:
            raise ValueError('The velocity fields should have the same shape.')
        # Compute perpendicular vectors.
        x_norm,y_norm,x_perp,y_perp,x_perp_all,y_perp_all=self.vertor_perp()
        # Define locations to interpolate interpolation.
        # !Important: 
        # x_perp,y_perp is defined in the center of the cells
        x = xr.DataArray(x_perp, dims=('transect','n'))
        y = xr.DataArray(y_perp, dims=('transect','n'))
        # Define limits to slice data.
        deltax    = 2*max((abs(x_perp[:,0]-x_perp[:,1])))
        slicevalx = [360+x_perp.min()-deltax,360+x_perp.max()+deltax]
        deltay    = 2*max((abs(y_perp[:,0]-y_perp[:,1])))
        slicevaly = [y_perp.min()-deltay,y_perp.max()+deltay]
        # Slice data to reduce memory issues.
        u = u.sel({'lon':slice(slicevalx[0],slicevalx[1]),'lat':slice(slicevaly[0],slicevaly[1])})
        v = v.sel({'lon':slice(slicevalx[0],slicevalx[1]),'lat':slice(slicevaly[0],slicevaly[1])})
        # Interpolate data using xarray, 
        # Note that fields can not contain nans
        # TO DO: Add support for data containing nans.
        print('Interpolating velocity fields')
        interp_u = u.interp(lon=360+x,lat=y).compute()
        interp_u = interp_u.where(interp_u!=0,np.nan)
        interp_v = v.interp(lon=360+x,lat=y).compute()
        interp_v = interp_v.where(interp_v!=0,np.nan)
        # Merge datasets
        self.interp_data=xr.merge([interp_u.to_dataset(name='u'), interp_v.to_dataset(name='v')])
        # Interpolate tracer fields to constrain transport.
        if tracerfile != None:
            print('Loadng and interpolating tracer')
            tracer = self.loaddata(file=tracerfile,var='Tracer',dataset=dataset,**xr_openmf_defaults)
            tracer = tracer.sel({'lon':slice(slicevalx[0],slicevalx[1]),'lat':slice(slicevaly[0],slicevaly[1])})
            interp_tracer = tracer.interp(lon=360+x,lat=y).compute()
            interp_tracer = interp_tracer.where(interp_tracer!=0,np.nan)
            self.interp_data   = xr.merge([interp_u.to_dataset(name='u'), interp_v.to_dataset(name='v'),
                                  interp_tracer.to_dataset(name='tracer')])
        # Save data.
        if save==True:
            self.interp_data.to_netcdf('./tmp_interp_transects.nc')
        return self.interp_data
    
    def depth_profiles(self,bottom_vel):
        '''
        
        '''
        # Maximum depth from interpolated field. 
        depth_index=self.interp_data.depth[np.isfinite(self.interp_data.u.where(abs(self.interp_data.u)>bottom_vel,np.nan).isel({'time':0})).argmin('depth')]
        # xr.DataArray 2 multiply with field.
        depth=(xr.zeros_like(self.interp_data.u.isel(time=0))+self.interp_data.depth)
        # Mask depth to only contain values larger than index.
        depth=depth.where(depth > depth_index,np.nan)
        # Delta depth to compute area
        delta_depth=depth.diff(dim='depth')
        return delta_depth
    
    def vel_magnitude(self):
        # Magnitude of interpolated vectors.
        magnitude = np.sqrt(self.interp_data.u**2+self.interp_data.v**2)
        return magnitude
    
    def dot_product(self):
        # Dot product between interpolated vectors and normal vector 
        # from perpendicular transect to the coast.
        return self.interp_data.u*self.x_norm[np.newaxis,np.newaxis,:,np.newaxis]+self.interp_data.v*self.y_norm[np.newaxis,np.newaxis,:,np.newaxis]
    
    def compute_transport(self,bottom_vel=1e-5):
        # Scalar projection of interpolated data
        dotproduct = self.dot_product()
        # Projected data over normal vectors to surface.
        u_normal   = dotproduct*self.x_norm[np.newaxis,np.newaxis,:,np.newaxis]
        v_normal   = dotproduct*self.y_norm[np.newaxis,np.newaxis,:,np.newaxis]
        # Area of each grid cell.
        dA = self.delta_area(bottom_vel)
        # Multiplication of vector sum and the dA. Flux integral.
        self.transport=(u_normal+v_normal)*dA
        return self.transport.sum(dim={'depth','n'})
    
    def delta_area(self,bottom_vel):
        # Compute perpendicular vectors.
        x_norm,y_norm,x_perp,y_perp,x_perp_all,y_perp_all=self.vertor_perp()
        # Depth at each section of the transect.
        delta_z=abs(self.depth_profiles(bottom_vel=bottom_vel))
        # Distance between lon,lat points of transect.
        delta_x=gsw.distance(x_perp_all,y_perp_all)
        return delta_z*delta_x
    
    def mask_transport(self,threshold,method='greater'):
        '''
        threshold [ float / list ]
            Threshold to scale transport with tracers used for tracer.
        method     [ string ] 
            'greater' will compute the transport for all the values larger 
                      than the threshold in the tracer field.
            'smaller' will compute the transport for all the values smaller 
                      than the threshold in the tracer field.
            'both' will compute the transport for all the values within 
                      the threshold interval in the tracer field.
        '''
        if type(threshold)==list:
            threshold=np.array(threshold)
        # TO DO: If u vertical grid != tracer vertical grid then interpolate tracer to velocity grid.
        if method=='smaller' and type(threshold)==float:
            scaled_transport=self.transport.where(self.interp_data.tracer.isel(depth=slice(0,-1))<threshold)
        elif method=='greater' and type(threshold)==float:
            scaled_transport=self.transport.where(self.interp_data.tracer.isel(depth=slice(0,-1))>threshold)
        elif method=='both' and type(threshold)==np.ndarray:
            scaled_transport=self.transport.where(self.interp_data.tracer.isel(depth=slice(0,-1))>threshold.min()).where(self.interp_data.tracer<threshold.max())
        else:
            raise ValueError('''Threshold must be an float or list/array in which the 
                             min and max value will define the threshold interval.''')
        return scaled_transport.sum(dim={'depth','n'})
    
    def loaddata(self,file=None,var='U',dataset=None,**kwargs):
        # Check if file or dataset is defined.
        if file == None and dataset==None:
            raise ValueError('''file should be the path to the netCDF files or 
                                dataset should contain a dataset with a variable 
                                containing the string defined as var.
                             ''')
        elif file != None and dataset == None:
            results = subprocess.check_output(['find', self.path, '-name', file])
            results=[s for s in results.decode('utf-8').split()]
            results.sort()
            data=xr.open_mfdataset(results,**kwargs)
        elif dataset != None and file == None:
            data=dataset
        else:
            raise ValueError('Only on of the arguments [file or dataset] can be defined.')
        # Extract variables from dataset
        varname= [key for key,items in data.data_vars.items()]
        # Rename variable for easier manipulation.
        if len(varname)==1:
            variable=data.rename({varname[0]:var})
        else:
            varname=[var for varn in varname if var in varn]
            variable=data.rename({varname[0]:var})
        # Extract only the variable of interest.
        data=variable[var]
        if type(data) != xr.core.dataarray.DataArray:
            raise ValueError('The provided data should be a xr.DataArray.')
        else:
            return data
        
    def vector_scale(self,index_perp,perp_angle):
        '''
        Scale vector to desired distance self.length
        '''
        # Scale perpendicular vector to distance self.length
        return np.squeeze((self.length*1000)/self.perpvecdist(index_perp,perp_angle))
    
    def vertor_perp(self,shift=0):
        # Nearest location of perpendicular vectors from coastline grid.
        index_perp=self.perploc()
        # Compute perpendicular vectors.
        perp_dict=self.perp2coast()
        # Scale perpendicular vector to desired distance self.length.
        scale=self.vector_scale(index_perp,perp_dict['Pvector']['angles'])
        # Gridded normal vector
        x_norm=(np.squeeze(np.linspace(0,scale,self.length//self.n)[:,np.newaxis]*self.x_norm)
                +self.coastline[index_perp][:,0]).T+shift
        y_norm=(np.squeeze(np.linspace(0,scale,self.length//self.n)[:,np.newaxis]*self.y_norm)
                +self.coastline[index_perp][:,1]).T
        # Gridded perpendicular vector at [x,y]
        x_perp_all=(np.squeeze(np.linspace(0,scale,self.length//self.n)[:,np.newaxis]*self.x_perp)
                +self.coastline[index_perp][:,0]).T+shift
        y_perp_all=(np.squeeze(np.linspace(0,scale,self.length//self.n )[:,np.newaxis]*self.y_perp)
                +self.coastline[index_perp][:,1]).T
        # Gridded perpendicular vector at [x+diff(x)/2,y+diff(y)/2]
        x_perp = x_perp_all[:,:-1]+np.diff(x_perp_all)/2
        y_perp = y_perp_all[:,:-1]+np.diff(y_perp_all)/2
        return x_norm,y_norm,x_perp,y_perp,x_perp_all,y_perp_all
        
    def plotperp_vect(self,shift=0,transect=None,**kwargs):
        '''
        transect [int] zooms in to the transect
        '''
        fig,ax = plt.subplots(1,1,figsize=(5,5),**kwargs)
        # Plot coastline
        plt.plot(self.coastline[:,0]+shift,self.coastline[:,1])
        # Compute perpendicular vectors.
        x_norm,y_norm,x_perp,y_perp,x_perp_all,y_perp_all=self.vertor_perp(shift)
        # Plot perpendicular vectors.
        plt.plot(x_norm.T,y_norm.T,'-r')
        plt.plot(x_perp.T,y_perp.T,'--k')
        # Zoom in into transect, usefule when the angle is significantly
        # different to n*{0,np.pi/2}.
        if transect != None:
            xdelta=2*abs(x_perp[transect][0]-x_perp[transect][1])
            plt.xlim(x_perp[transect].min()-xdelta,x_perp[transect].max()+xdelta)
            ydelta=2*abs(y_perp[transect][0]-y_perp[transect][1])
            plt.ylim(y_perp[transect].min()-ydelta,y_perp[transect].max()+ydelta)
        plt.gca().set_aspect('equal', adjustable='box')
        return fig,ax