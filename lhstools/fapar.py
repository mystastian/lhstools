from .benchmark import Benchmark
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import numpy as np


class FAPAR(Benchmark):
    """FAPAR benchmark"""
    def __init__(self,config_name='config.ini',init=False,*args,**kwargs):
        #Shared parameters as attributes
        super(FAPAR,self).__init__(config_name,init,*args,**kwargs)
        #FAPAR parameters,
        super(FAPAR,self).import_config(config_name,'FAPAR')
        #Load observation file and rename/reorder dimensions
        self.obs=xr.open_dataset(self.fnobs,decode_times=False).FAPAR
        self.obs=self.obs.rename({'TAX':'TIME'})
        self.obs=self.obs.transpose('TIME','LATITUDE','LONGITUDE')
        #Overwrite timeaxis with datetime object
        self.timeaxis=pd.date_range('1997-09-01',freq='m',periods=106)
        self.obs['TIME']=self.timeaxis
        self.obs_mean=self.obs.mean(dim='TIME')
        #Monthly mean
        self.obs_monthly=self.obs.groupby('TIME.month').mean(dim='TIME')
        #Load par from control run
        self.par=xr.open_dataset(self.fnctrl_m).par.isel(grid_only=0).drop('grid_only')
        self.tconstraint=(self.par.TIME>1997.65) * (self.par.TIME<2006.5)
        self.par=self.par[self.tconstraint]
        self.par['TIME']=self.timeaxis
        self.name='fapar'
        #Load area from control run, stack it for seasonal mean
        self.area=xr.open_dataset(self.fnctrl).area
        self.area_monthly=xr.concat([self.area]*12,dim='month')
        self.area_monthly.month.values=range(1,13)

    def calc_stats(self,memberid):
        "Returns stats (error and variance) of a member id"
        #NEW: error of every month
        fapar_sim=self.get_sim(memberid).groupby('TIME.month').mean(dim='TIME')
        return super(FAPAR,self).calc_metric(fapar_sim,self.obs_monthly,self.area_monthly)
        #OLD: Simple mean
        # fapar_sim=self.get_sim(memberid).mean(axis=0)
        # #Error according to specified metric
        # return Benchmark.calc_metric(self,fapar_sim,self.obs_mean,self.area)


    def get_sim(self,memberid):
        "returns mean simulated fapar"
        fnmember=self.path2cdf+memberid+'_m.cdf'
        #Same time range as in SeaWifs dataset
        #Get apar and luarea
        apar=xr.open_dataset(fnmember).apar[self.tconstraint]
        apar['TIME']=self.timeaxis
        luarea=xr.open_dataset(fnmember).lu_area[self.tconstraint]
        luarea['TIME']=self.timeaxis
        #Calculate simulated FAPAR
        return (apar*luarea).sum(dim='landuse')/luarea.sum(dim='landuse')/self.par

    def plot_member(self,memberid):
        "Plots given member and observed fapar"
        cm='bwr'
        projection=ccrs.PlateCarree()
        #prepare data
        #side by side
        fapar_sim=self.get_sim(memberid)
        diff=self.obs-fapar_sim
        fig1=plt.figure()
        ax = plt.subplot(211,projection=projection)
        ax.coastlines()
        self.obs_mean.plot.pcolormesh(ax=ax,vmin=0,vmax=0.7,cmap=cm)
        ax.set_title('SeaWiFS FAPAR on 1x1 grid, mean(09.1997-06.2006)')
        ax = plt.subplot(212,projection=ccrs.PlateCarree())
        ax.coastlines()
        fapar_sim.mean(dim='TIME').plot.pcolormesh(ax=ax,vmin=0,vmax=0.7,cmap=cm)
        ax.set_title('Control LPJ FAPAR, mean(09.1997-06.2006)')
        #Difference
        fig2=plt.figure()
        ax = plt.axes(projection=projection)
        ax.coastlines()
        diff.mean(dim='TIME').plot.pcolormesh(ax=ax)
        #Difference every month
        fig3,axs = plt.subplots(4,3,figsize=(16,16),subplot_kw={'projection':projection})

        ax_cbar=fig3.add_axes([0.92, 0.1, 0.02, 0.2])
        fig3.subplots_adjust(right=0.9)
        axs=np.ndarray.flatten(axs)
        for m in range(12):
            axs[m].coastlines()
            im=diff.groupby('TIME.month').mean(dim='TIME').sel(month=m+1).plot.pcolormesh(ax=axs[m],vmin=-0.5,vmax=0.5,add_colorbar=False,cmap=cm)
        plt.colorbar(im,cax=ax_cbar)

        #Zonal means
        fig3,ax=plt.subplots()
        self.obs_mean.mean(dim='LONGITUDE').plot(label='SeaWiFS')
        fapar_sim.mean(dim=['TIME','LONGITUDE']).plot(label='Simulation: {}'.format(memberid))
        ax.set_title('Zonal and temporal FAPAR mean')
        plt.legend(fancybox=True)
        self.obs_mean
        return fig1,fig2,fig3



