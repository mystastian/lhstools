from .benchmark import Benchmark
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import numpy as np
from lhstools.utils import discrete_cmap

class BiomassKeith(Benchmark):
    """Biomass Keith et. al. 2009 benchmark"""
    def __init__(self,config_name='config.ini',init=False,*args,**kwargs):
        #Shared parameters as attributes
        Benchmark.__init__(self,config_name,init,*args,**kwargs)
        #FAPAR parameters,
        Benchmark.import_config(self,config_name,'BiomassKeith')
        #Load observation file and look for corresponding lat/lon coordinates in control file
        self.obs=pd.read_csv(self.fnobs,delim_whitespace=True,skiprows=11)
        #Column 'abg' (Aboveground biomass) Conversion of 100 for tC/ha -> gC/m^2
        self.obs['vegc']=self.obs['agb']*100
        ctrl=xr.open_dataset(self.fnctrl)
        self.obs['sim_lat']=[float(ctrl.LATITUDE.sel(LATITUDE=lat,method='nearest')) for lat in self.obs['lat']]
        self.obs['sim_lon']=[float(ctrl.LONGITUDE.sel(LONGITUDE=lon,method='nearest')) for lon in self.obs['lon']]
        #Average observation that are in the same model grid cell.
        self.obs_ave=self.obs.groupby(['sim_lat','sim_lon']).mean()
        self.obs_ave_num=self.obs.groupby(['sim_lat','sim_lon']).count()
        self.name='biomasskeith'

    def calc_stats(self,memberid):
        "Returns stats (error and variance) of a member id"
        #NEW: error of every month
        vegc,dump=self.get_sim(memberid)
        return Benchmark.calc_metric(self,vegc,self.obs_ave['vegc'],weight=None)


    def get_sim(self,memberid):
        "Returns vegc of a member at the observation cells as a pandas series and the averaged cdf"
        fnmember=self.path2cdf+memberid+'.cdf'
        vegc_cdf=xr.open_dataset(fnmember).vegcarbon
        #Sum over Nat and Secd LU Class if gross LU
        if 'Secondary' in self.lunames:
            luarea_cdf=xr.open_dataset(fnmember).lu_area
            vegc_cdf=vegc_cdf.sel(landuse=[1,2])
            luarea=luarea_cdf.sel(landuse=[1,2])
            #Not sure if actually correct this way..
            vegc_cdf=(vegc_cdf*luarea).sum(dim='landuse')/luarea.sum(dim='landuse')
        else:
            vegc_cdf=vegc_cdf.sel(landuse=1)
        #Time constraint [Rather arbitrary right now..]
        tconstraint=(vegc_cdf.TIME>1950.)*(vegc_cdf.TIME<2000.)
        vegc_cdf=vegc_cdf[tconstraint].mean(dim='TIME')
        #Panda series with obs_ave index
        vegc=pd.Series(index=self.obs_ave.index)
        for coord,dump in vegc.iteritems():
            vegc.ix[coord]=float(vegc_cdf.sel(LATITUDE=coord[0],LONGITUDE=coord[1]))
        return vegc,vegc_cdf



    def plot_member(self,memberid):
        "Plots given member and observed Keith Biomass"
        #Plot setings
        cm='jet'
        levels=11
        dcm=discrete_cmap(cm,levels)
        maxvegc=40000.

        #Import
        vegc_sites,vegc_map=self.get_sim(memberid)

        #Plot map and extract simulation values
        fig1,ax = plt.subplots(figsize=(10,4),subplot_kw={'projection':ccrs.PlateCarree()})
        CS=vegc_map.plot.contourf(cmap=cm,vmax=maxvegc,levels=levels,add_colorbar=False)

        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(r'Vegetation carbon [gC m$^2$ yr$^{-1}$]')
        for coord,svegc in self.obs['vegc'].iteritems():
            color=dcm(int(np.floor((svegc/maxvegc*levels))))
            ax.scatter(self.obs['lon'].ix[coord],self.obs['lat'].ix[coord],200,color=color,marker='.',edgecolor='k')
            if svegc<1.0:
                ax.scatter(lon,lat,200,color=color,marker='*',edgecolor='k')
        ax.set_title('Simulated Vegetation Carbon ('+memberid+') vs BiomassKeith Class A measurements')
        ax.coastlines()

        #Scatter Plot
        fig2,ax = plt.subplots()
        plt.scatter(x=vegc_sites,y=self.obs_ave['vegc'])
        plt.scatter(x=vegc_sites,y=self.obs_ave['vegc'])
        plt.scatter(x=vegc_sites,y=self.obs_ave['vegc'])

        ax.set_title('Simulated Vegetation Carbon ('+memberid+') vs BiomassKeith Class A measurements')
        ax.set_xlabel(r'Vegetation Carbon simulated [gC m$^2$ yr$^{-1}$]')
        ax.set_ylabel(r'Vegetation Carbon observed (BiomassKeith) [gC m$^2$ yr$^{-1}$]')
        ax.plot([0,maxvegc],[0,maxvegc],color='k')
        ax.set_xlim(-5,maxvegc)
        ax.set_ylim(-5,maxvegc)

        text='Correlation (R^2): {0}'.format(vegc_sites.corr(self.obs_ave['vegc']))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
        return fig1,fig2


