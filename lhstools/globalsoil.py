from .benchmark import Benchmark
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import numpy as np
from lhstools.utils import discrete_cmap

class GlobalSoil(Benchmark):
    """Global Soil carbon inventory"""
    def __init__(self,config_name='config.ini',init=False,*args,**kwargs):
        #Shared parameters as attributes
        Benchmark.__init__(self,config_name,init,*args,**kwargs)
        #FAPAR parameters,
        Benchmark.import_config(self,config_name,'GlobalSoil')
        #Get observations
        self.get_obs()
        #name
        self.name='GlobalSoil'

    def get_obs(self):
        """set data according to
        IPCC Working Group 1
        """
        self.obs=pd.Series(index=['SoilCarbon'])

        #From IPCC Infog
        self.obs['SoilCarbon']=1950
        self.sigma_obs=450



    def calc_stats(self,memberid):
        "Returns stats (error and variance) of a member id"
        #NEW: error of every month
        soilc=self.get_sim(memberid)
        return Benchmark.calc_metric(self,soilc,self.obs,weight=None,sigma_obs=self.sigma_obs)


    def get_sim(self,memberid):
        "Returns the  land emissions from the totc ascii of a member as a pandas dataframe"
        fnmember=self.path2ascii+'trans_'+memberid+'.soilc.out'
        soilc=Benchmark.read_ascii(self,fnmember)['Total']
        #average over whole simulated period...
        meansoil=pd.Series(index=['SoilCarbon'])
        meansoil['SoilCarbon']=soilc.mean(axis=0)
        return meansoil

    def plot_hist(self):
        "Pltos a histogram of all the vegetation for all members"
        #Import all the data
        soil=pd.DataFrame(index=self.obs.index)
        for member in self.members:
            soil[member]=self.get_sim(member)
        # ax=uptake.transpose().hist()
        fig,ax=plt.subplots(figsize=(12,8))
        ax=soil.transpose().hist(ax=ax, fc='lightblue', histtype='stepfilled', alpha=0.3, normed=True,label='Norm. Histogram')
        plt.suptitle('Soil Carbon')
        #Add observations and kde
        soil.transpose()['SoilCarbon'].plot(kind='kde',label='KDE',ax=ax[0],bw_method=0.5,color='lightblue')
        ax[0].plot([self.obs['SoilCarbon'],self.obs['SoilCarbon']],ax[0].get_ylim(),color='red',linewidth=2,label='IPCC')
        ax[0].axvspan(self.obs['SoilCarbon']-self.sigma_obs,self.obs['SoilCarbon']+self.sigma_obs,color='red',alpha=0.3,label=r'IPCC $\pm \sigma$')
        ax[0].legend(loc='upper left',fancybox=True,framealpha=0.8)
        return fig,ax










