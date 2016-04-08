from benchmark import Benchmark
import xarray as xr
import cartopy.crs as ccrs

class FAPAR(Benchmark):
    """FAPAR benchmark"""
    def __init__(self,config_name='config.ini'):
        #Shared parameters as attributes
        Benchmark.__init__(self,config_name)
        #FAPAR parameters,
        Benchmark.import_config(self,config_name,'FAPAR')
        #Load observation file
        self.obs=xr.open_dataset(self.fnobs,decode_times=False).FAPAR
        self.obs_mean=self.obs.mean(axis=0).transpose()
        #Load par from seperate file
        self.par=xr.open_dataset(self.fnctrl_m).par.isel(grid_only=0)

    def calc_error(self,memberid):
        "Returns error of a member id"
        fnmember=self.path2cdf+memberid+'_m.cdf'
        #Same time range as in SeaWifs dataset
        tconstraint=(self.par.TIME>1997.65) * (self.par.TIME<2006.5)
        #Get apar and luarea
        apar=xr.open_dataset(fnmember).apar[tconstraint]
        luarea=xr.open_dataset(fnmember).lu_area[tconstraint]
        #Calculate simulated FAPAR
        fapar_sim=(apar*luarea).sum(axis=1)/self.par[tconstraint]
        #Error according to specified metric
        return Benchmark.calc_metric(self,fapar_sim.mean(axis=0),self.obs_mean)




