from .benchmark import Benchmark
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import numpy as np


class ETMap(Benchmark):
    """Evapotranspiration Map benchmark from LandFlux-Eval"""
    def __init__(self,config_name='config.ini',init=False,*args,**kwargs):
        #Shared parameters as attributes
        super(ETMap,self).__init__(config_name,init,*args,**kwargs)
        #FAPAR parameters,
        super(ETMap,self).import_config(config_name,'ETMap')
        #Load observation file and rename/reorder dimensions
        self.obs=xr.open_dataset(self.fnobs).ET_mean*365
        self.obs=self.obs.rename({'lon':'LONGITUDE','lat':'LATITUDE'})

        #Load area from control run, stack it for seasonal mean
        self.area=xr.open_dataset(self.fnctrl).area
        #Name
        self.name='ETMap'

    def calc_stats(self,memberid):
        "Returns stats (error and variance) of a member id"
        etransp_sim=self.get_sim(memberid)
        return super(ETMap,self).calc_metric(etransp_sim,self.obs,self.area)


    def get_sim(self,memberid):
        "Get mean evapotranspiration from 1989-2005"
        fnmember=self.path2cdf+memberid+'.cdf'
        #Open Datasets
        evap=xr.open_dataset(fnmember).evap
        transp=xr.open_dataset(fnmember).transp
        etransp=evap+transp
        luarea=xr.open_dataset(fnmember).lu_area
        #Time constraint (+0.6 to count endyear as well..)
        tconstraint=(evap.TIME>float(self.startyear)) * (evap.TIME<float(self.endyear)+0.6)
        sum_etransp=(etransp[tconstraint]*luarea[tconstraint]).sum(dim='landuse',skipna=False).mean(dim='TIME',skipna=False)
        return sum_etransp


    def plot_member(self,memberid):
        "Plots given member and observed evapotranspiration"
        cm='bwr'
        projection=ccrs.PlateCarree()
        #prepare data
        #side by side
        sim=self.get_sim(memberid)
        diff=self.obs-sim
        fig1=plt.figure()
        ax = plt.subplot(211,projection=projection)
        ax.coastlines()
        CS=self.obs.plot.contourf(ax=ax,vmin=0,vmax=2000,cmap='Reds',add_colorbar=False)
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(r'Evapotranspiration [mm yr^{-1]]')
        ax.set_title('LandFluxEVAL evapotranspiration')
        ax = plt.subplot(212,projection=ccrs.PlateCarree())
        ax.coastlines()
        CS=sim.plot.contourf(ax=ax,vmin=0,vmax=2000,cmap='Reds',add_colorbar=False)
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(r'Evapotranspiration [mm yr^{-1]]')
        ax.set_title('LPX '+memberid+' evapotranspiration, mean(+'+self.startyear+','+self.endyear+')')
        plt.tight_layout()
        #Difference
        fig2=plt.figure(figsize=(8,4))
        ax=plt.subplot(111,projection=projection)
        ax.set_title('LandFluxEVAL evapotranspiration - LPX '+memberid+' evapotranspiration (difference)')
        cbar.ax.set_ylabel(r'Evapotranspiration [mm yr^{-1]]')
        ax.coastlines()
        diff.plot.contourf(ax=ax,levels=18)
        return fig1,fig2



