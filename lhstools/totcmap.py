from .benchmark import Benchmark
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import numpy as np


class TotCMap(Benchmark):
    """TotalCarbon Map benchmark from Carvalhais et. al. (2014)"""
    def __init__(self,config_name='config.ini',init=False,*args,**kwargs):
        #Shared parameters as attributes
        super(TotCMap,self).__init__(config_name,init,*args,**kwargs)
        #FAPAR parameters,
        super(TotCMap,self).import_config(config_name,'TotCMap')
        #Load observation file and rename/reorder dimensions
        self.obs=xr.open_dataset(self.fnobs).cTotal
        #Load area from control run, stack it for seasonal mean
        self.area=xr.open_dataset(self.fnctrl).area
        #Name
        self.name='TotCMap'

    def calc_stats(self,memberid):
        "Returns stats (error and variance) of a member id"
        totc_sim=self.get_sim(memberid)
        return super(TotCMap,self).calc_metric(totc_sim,self.obs,self.area)


    def get_sim(self,memberid):
        "returns mean simulated vegc+soilc+litterc+exuc (No prodc) between 1982 and 2005 (As specified in Carvalhais, last page)"
        fnmember=self.path2cdf+memberid+'.cdf'
        #Open Datasets
        soil=xr.open_dataset(fnmember).soilcarbon
        litter_ag=xr.open_dataset(fnmember).littercarbon_ag
        litter_bg=xr.open_dataset(fnmember).littercarbon_bg
        vegc=xr.open_dataset(fnmember).vegcarbon
        exuc=xr.open_dataset(fnmember).exucarbon
        luarea=xr.open_dataset(fnmember).lu_area
        #Time constraint
        tconstraint=(soil.TIME>1982.1) * (soil.TIME<2005.9)
        sum_soil=(soil[tconstraint]*luarea[tconstraint]).sum(dim='landuse',skipna=False).mean(dim='TIME',skipna=False)
        sum_litter_ag=(litter_ag[tconstraint]*luarea[tconstraint]).sum(dim='landuse',skipna=False).mean(dim='TIME',skipna=False)
        sum_litter_bg=(litter_bg[tconstraint]*luarea[tconstraint]).sum(dim='landuse',skipna=False).mean(dim='TIME',skipna=False)
        sum_vegc=(vegc[tconstraint]*luarea[tconstraint]).sum(dim='landuse',skipna=False).mean(dim='TIME',skipna=False)
        sum_exuc=(exuc[tconstraint]*luarea[tconstraint]).sum(dim='landuse',skipna=False).mean(dim='TIME',skipna=False)
        return sum_soil+sum_litter_ag+sum_litter_bg+sum_litter_bg+sum_vegc+sum_exuc


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
        cbar.ax.set_ylabel(r'Totalcarbon [gC m^{-2]]')
        ax.set_title('Carvalhais et. al Totalcarbon on 1x1 grid')
        ax = plt.subplot(212,projection=ccrs.PlateCarree())
        ax.coastlines()
        CS=sim.plot.contourf(ax=ax,vmin=0,vmax=100000,cmap='Reds',add_colorbar=False)
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(r'Totalcarbon [gC m^{-2]]')
        ax.set_title('LPX '+memberid+' Totalcarbon, mean(1982,2005)')
        plt.tight_layout()
        #Difference
        fig2=plt.figure(figsize=(8,4))
        ax=plt.subplot(111,projection=projection)
        ax.set_title('Carvalhais Totalcarbon - LPX '+memberid+' Totalcarbon difference')
        ax.coastlines()
        diff.plot.contourf(ax=ax,levels=18)
        return fig1,fig2



