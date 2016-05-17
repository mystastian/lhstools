from .benchmark import Benchmark
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import numpy as np
from lhstools.utils import discrete_cmap

class FLUXNET(Benchmark):
    """FLUXNET NPP benchmark"""
    def __init__(self,config_name='config.ini',init=False,*args,**kwargs):
        #Shared parameters as attributes
        Benchmark.__init__(self,config_name,init,*args,**kwargs)
        #FAPAR parameters,
        Benchmark.import_config(self,config_name,'FLUXNET')
        #Load observation file and look for corresponding lat/lon coordinates in control file
        self.obs=pd.read_csv(self.fnobs,delim_whitespace=True,header=4,na_values=-9999)
        ctrl=xr.open_dataset(self.fnctrl)
        #Use tnpp2, where available, otherwise tnpp1, if none is available drop the column
        self.obs['npp']=self.obs['tnpp2'].fillna(self.obs['tnpp1'])
        self.obs=self.obs[np.isfinite(self.obs['npp'])]
        self.obs.index=range(len(self.obs))
        self.obs['lat_sim']=[float(ctrl.LATITUDE.sel(LATITUDE=lat,method='nearest')) for lat in self.obs['lat']]
        self.obs['lon_sim']=[float(ctrl.LONGITUDE.sel(LONGITUDE=lon,method='nearest')) for lon in self.obs['lon']]
        #Average observation that are in the same model grid cell.
        self.obs_ave=self.obs.groupby(['lat_sim','lon_sim']).mean()
        self.obs_ave_num=self.obs.groupby(['lat_sim','lon_sim']).count()
        self.name='fluxnet'

    def calc_stats(self,memberid):
        "Returns stats (error and variance) of a member id"
        #NEW: error of every month
        npp,dump=self.get_sim(memberid)
        return Benchmark.calc_metric(self,npp,self.obs_ave['npp'],weight=None)


    def get_sim(self,memberid):
        "Returns npp of a member at the observation cells as a pandas series and the averaged cdf"
        fnmember=self.path2cdf+memberid+'.cdf'
        npp_cdf=xr.open_dataset(fnmember).npp
        luarea_cdf=xr.open_dataset(fnmember).lu_area
        #Sum over Nat and Secd LU Class
        npp_cdf=npp_cdf.sel(landuse=[1,2])
        luarea=luarea_cdf.sel(landuse=[1,2])
        npp_cdf=(npp_cdf*luarea).sum(dim='landuse')/luarea.sum(dim='landuse')
        #Time constraint
        tconstraint=(npp_cdf.TIME>1931.)*(npp_cdf.TIME<1997.)
        npp_cdf=npp_cdf[tconstraint].mean(dim='TIME')
        #Panda series with obs_ave index
        npp=pd.Series(index=self.obs_ave.index)
        for coord,dump in npp.iteritems():
            npp.ix[coord]=float(npp_cdf.sel(LATITUDE=coord[0],LONGITUDE=coord[1]))
        return npp,npp_cdf



    def plot_member(self,memberid):
        "Plots given member and observed fapar"
        #Plot setings
        cm='jet'
        levels=11
        dcm=discrete_cmap(cm,levels)
        maxNPP=1500.

        #Import
        npp_sites,npp_map=self.get_sim(memberid)

        #Plot map and extract simulation values
        fig1,ax = plt.subplots(figsize=(10,4),subplot_kw={'projection':ccrs.PlateCarree()})
        npp_map.plot.contourf(cmap=cm,vmax=maxNPP,levels=levels)
        for coord,sNPP in self.obs['npp'].iteritems():
            npp_site=npp_sites[self.obs['lat_sim'][coord],self.obs['lon_sim'][coord]]
            color=dcm(int(np.floor((sNPP/maxNPP*levels))))
            if (npp_site<0.1) or (np.isnan(npp_site)):
                ax.scatter(self.obs['lon'].ix[coord],self.obs['lat'].ix[coord],200,color=color,marker='*',edgecolor='k')
            else:
                ax.scatter(self.obs['lon'].ix[coord],self.obs['lat'].ix[coord],200,color=color,marker='.',edgecolor='k')
                # pass
        ax.set_title('Simulated NPP ('+memberid+') vs FLUXNET Class A measurements')
        ax.coastlines()

        #Scatter Plot
        fig2,ax = plt.subplots()
        plt.scatter(x=npp_sites,y=self.obs_ave['npp'])
        ax.set_title('Simulated NPP ('+memberid+') vs FLUXNET Class A measurements')
        ax.set_xlabel(r'NPP simulated [gc/m^2/yr]')
        ax.set_ylabel(r'NPP observed (FLUXNET) [gc/m^2/yr]')
        ax.plot([0,maxNPP],[0,maxNPP],color='k')
        ax.set_xlim(-5,maxNPP)
        ax.set_ylim(-5,maxNPP)
        #Fit (Not currently used)
        # model=pd.ols(x=obs.sim_npp,y=obs.npp)
        # x=np.linspace(0,maxNPP,76)
        # p=np.poly1d(model.beta)
        # plt.plot(x,p(x),label='Fit: y={0:.3f}x+{1:.3f}\n R$^2$={2:.3f}'.format(model.beta.x,model.beta.intercept,model.r2))


        text='Correlation (R^2): {0}'.format(npp_sites.corr(self.obs_ave['npp']))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
        return fig1,fig2


