from .benchmark import Benchmark
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import numpy as np
from lhstools.utils import discrete_cmap

class GlobalVeg(Benchmark):
    """Global vegetation carbon inventory"""
    def __init__(self,config_name='config.ini',init=False,*args,**kwargs):
        #Shared parameters as attributes
        Benchmark.__init__(self,config_name,init,*args,**kwargs)
        #FAPAR parameters,
        Benchmark.import_config(self,config_name,'GlobalVeg')
        #Get observations
        self.get_obs()
        #name
        self.name='GlobalVeg'

    def get_obs(self):
        """set data according to
        IPCC Working Group 1
        """
        self.obs=pd.Series(index=['VegCarbon'])

        #From IPCC Infog
        self.obs['VegCarbon']=550
        self.sigma_obs=100



    def calc_stats(self,memberid):
        "Returns stats (error and variance) of a member id"
        #NEW: error of every month
        vegc=self.get_sim(memberid)
        return Benchmark.calc_metric(self,vegc,self.obs,weight=None,sigma_obs=self.sigma_obs)


    def get_sim(self,memberid):
        "Returns the  land emissions from the totc ascii of a member as a pandas dataframe"
        fnmember=self.path2ascii+'trans_'+memberid+'.vegc.out'
        vegc=Benchmark.read_ascii(self,fnmember)['Total']
        #average over whole simulated period...
        meanveg=pd.Series(index=['VegCarbon'])
        meanveg['VegCarbon']=vegc.mean(axis=0)
        return meanveg

    def plot_hist(self):
        "Pltos a histogram of all the vegetation for all members"
        #Import all the data
        veg=pd.DataFrame(index=self.obs.index)
        for member in self.members:
            veg[member]=self.get_sim(member)
        # ax=uptake.transpose().hist()
        fig,ax=plt.subplots(figsize=(12,8))
        ax=veg.transpose().hist(ax=ax, fc='lightblue', histtype='stepfilled', alpha=0.3, normed=True,label='Norm. Histogram')
        plt.suptitle('Vegetation Carbon')
        #Add observations and kde
        veg.transpose()['VegCarbon'].plot(kind='kde',label='KDE',ax=ax[0],bw_method=0.5,color='lightblue')
        ax[0].plot([self.obs['VegCarbon'],self.obs['VegCarbon']],ax[0].get_ylim(),color='red',linewidth=2,label='IPCC')
        ax[0].axvspan(self.obs['VegCarbon']-self.sigma_obs,self.obs['VegCarbon']+self.sigma_obs,color='red',alpha=0.3,label=r'IPCC $\pm \sigma$')

        ax[0].legend(loc='upper left',fancybox=True,framealpha=0.8)

        return fig,ax










