from .benchmark import Benchmark
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import numpy as np
from lhstools.utils import discrete_cmap

class EMDI(Benchmark):
    """EMDI NPP benchmark"""
    def __init__(self,config_name='config.ini',init=False,*args,**kwargs):
        #Shared parameters as attributes
        Benchmark.__init__(self,config_name,init,*args,**kwargs)
        #FAPAR parameters,
        Benchmark.import_config(self,config_name,'EMDI')
        #Load observation file and look for corresponding lat/lon coordinates in control file
        self.obs=pd.read_csv(self.fnobs)
        ctrl=xr.open_dataset(self.fnctrl)
        self.obs['lat']=[float(ctrl.LATITUDE.sel(LATITUDE=lat,method='nearest')) for lat in self.obs['LAT_DD']]
        self.obs['lon']=[float(ctrl.LONGITUDE.sel(LONGITUDE=lon,method='nearest')) for lon in self.obs['LONG_DD']]
        #Average observation that are in the same model grid cell.
        self.obs_ave=self.obs.groupby(['lat','lon']).mean()
        self.obs_ave_num=self.obs.groupby(['lat','lon']).count()
        self.name='emdi'

    def calc_stats(self,memberid):
        "Returns stats (error and variance) of a member id"
        #NEW: error of every month
        npp,dump=self.get_sim(memberid)
        return Benchmark.calc_metric(self,npp,self.obs_ave['TNPP_C'],weight=None)


    def get_sim(self,memberid):
        "Returns npp of a member at the observation cells as a pandas series and the averaged cdf"
        fnmember=self.path2cdf+memberid+'.cdf'
        npp_cdf=xr.open_dataset(fnmember).npp
        luarea_cdf=xr.open_dataset(fnmember).lu_area
        #Sum over Nat and Secd LU Class if gross LU
        if 'Secondary' in self.lunames:
            npp_cdf=npp_cdf.sel(landuse=[1,2])
            luarea=luarea_cdf.sel(landuse=[1,2])
            #Not sure if actually correct this way..
            npp_cdf=(npp_cdf*luarea).sum(dim='landuse')/luarea.sum(dim='landuse')
        else:
            npp_cdf=npp_cdf.sel(landuse=1)
        #Time constraint
        tconstraint=(npp_cdf.TIME>1931.)*(npp_cdf.TIME<1997.)
        npp_cdf=npp_cdf[tconstraint].mean(dim='TIME')
        #Panda series with obs_ave index
        npp=pd.Series(index=self.obs_ave.index)
        for coord,dump in npp.iteritems():
            npp.ix[coord]=float(npp_cdf.sel(LATITUDE=coord[0],LONGITUDE=coord[1]))
        return npp,npp_cdf



    def plot_member(self,memberid):
        "Plots given member and observed EMDI NPP"
        #Plot setings
        cm='jet'
        levels=11
        dcm=discrete_cmap(cm,levels)
        maxNPP=1500.

        #Import
        npp_sites,npp_map=self.get_sim(memberid)

        #Plot map and extract simulation values
        fig1,ax = plt.subplots(figsize=(10,4),subplot_kw={'projection':ccrs.PlateCarree()})
        CS=npp_map.plot.contourf(cmap=cm,vmax=maxNPP,levels=levels,add_colorbar=False)

        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(r'NPP [gC m$^2$ yr$^{-1}$]')
        for coord,sNPP in self.obs['TNPP_C'].iteritems():
            color=dcm(int(np.floor((sNPP/maxNPP*levels))))
            ax.scatter(self.obs['LONG_DD'].ix[coord],self.obs['LAT_DD'].ix[coord],200,color=color,marker='.',edgecolor='k')
            if sNPP<1.0:
                ax.scatter(lon,lat,200,color=color,marker='*',edgecolor='k')
        ax.set_title('Simulated NPP ('+memberid+') vs EMDI Class A measurements')
        ax.coastlines()

        #Scatter Plot
        fig2,ax = plt.subplots()
        plt.scatter(x=npp_sites,y=self.obs_ave['TNPP_C'])
        plt.scatter(x=npp_sites,y=self.obs_ave['TNPP_C'])
        plt.scatter(x=npp_sites,y=self.obs_ave['TNPP_C'])

        ax.set_title('Simulated NPP ('+memberid+') vs EMDI Class A measurements')
        ax.set_xlabel(r'NPP simulated [gC m$^2$ yr$^{-1}$]')
        ax.set_ylabel(r'NPP observed (EMDI) [gC m$^2$ yr$^{-1}$]')
        ax.plot([0,maxNPP],[0,maxNPP],color='k')
        ax.set_xlim(-5,maxNPP)
        ax.set_ylim(-5,maxNPP)
        #Fit (Not currently used)
        # model=pd.ols(x=obs.sim_npp,y=obs.TNPP_C)
        # x=np.linspace(0,maxNPP,76)
        # p=np.poly1d(model.beta)
        # plt.plot(x,p(x),label='Fit: y={0:.3f}x+{1:.3f}\n R$^2$={2:.3f}'.format(model.beta.x,model.beta.intercept,model.r2))


        text='Correlation (R^2): {0}'.format(npp_sites.corr(self.obs_ave['TNPP_C']))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
        return fig1,fig2


