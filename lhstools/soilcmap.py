from .benchmark import Benchmark
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import numpy as np


class SoilCMap(Benchmark):
    """SoilCarbon Map benchmark from Carvalhais et. al. (2014)"""
    def __init__(self,config_name='config.ini',init=False,*args,**kwargs):
        #Shared parameters as attributes
        super(SoilCMap,self).__init__(config_name,init,*args,**kwargs)
        #FAPAR parameters,
        super(SoilCMap,self).import_config(config_name,'SoilCMap')
        #Load observation file and rename/reorder dimensions
        self.obs=xr.open_dataset(self.fnobs).cSoilTotal
        #Load area from control run, stack it for seasonal mean
        self.area=xr.open_dataset(self.fnctrl).area
        #Name
        self.name='SoilCMap'

    def calc_stats(self,memberid):
        "Returns stats (error and variance) of a member id"
        soilc_sim=self.get_sim(memberid)
        return super(SoilCMap,self).calc_metric(soilc_sim,self.obs,self.area)


    def get_sim(self,memberid):
        "returns mean simulated soilcarbon between 1982 and 2005 (As specified in Carvalhais, last page)"
        fnmember=self.path2cdf+memberid+'.cdf'
        #Open Dataset
        sim=xr.open_dataset(fnmember).soilcarbon
        luarea=xr.open_dataset(fnmember).lu_area
        #Time constraint
        tconstraint=(sim.TIME>1982.1) * (sim.TIME<2005.9)
        sum_sim=(sim[tconstraint]*luarea[tconstraint]).sum(dim='landuse',skipna=False).mean(dim='TIME',skipna=False)
        return sum_sim


    def plot_member(self,memberid):
        "Plots given member and observed fapar"
        cm='bwr'
        projection=ccrs.PlateCarree()
        #prepare data
        #side by side
        sim=self.get_sim(memberid)
        diff=self.obs-sim
        fig1=plt.figure()
        ax = plt.subplot(211,projection=projection)
        ax.coastlines()
        CS=self.obs.plot.contourf(ax=ax,vmin=0,vmax=100000,cmap='Reds',add_colorbar=False)
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(r'Soilcarbon [gC m^{-2]]')
        ax.set_title('Carvalhais et. al soilcarbon on 1x1 grid')
        ax = plt.subplot(212,projection=ccrs.PlateCarree())
        ax.coastlines()
        CS=sim.plot.contourf(ax=ax,vmin=0,vmax=100000,cmap='Reds',add_colorbar=False)
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(r'Soilcarbon [gC m^{-2]]')
        ax.set_title('LPX '+memberid+' soilcarbon, mean(1982,2005)')
        plt.tight_layout()
        #Difference
        fig2=plt.figure()
        ax = plt.axes(projection=projection)
        ax.set_title('Carvalhais Soilcarbon - LPX '+memberid+' Soilcarbon difference')
        ax.coastlines()
        diff.plot.contourf(ax=ax,levels=9)
        return fig1,fig2



