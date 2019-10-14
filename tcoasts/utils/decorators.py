import functools
import warnings
import os
import xarray as xr

def _file_exists(func):
    @functools.wraps(func)
    def wrap(self,*args, **kwargs):
        if os.path.isfile(self.tmpfile):
            warnings.filterwarnings('default',module='tcoasts')
            warnings.warn('Loading previous saved data.', Warning)
            interp=xr.open_dataset('./tmp_interp_transects.nc')
            self.interp_data=interp.load()
        else:
            func(self, *args, **kwargs)
    return wrap