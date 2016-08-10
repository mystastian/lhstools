from .benchmark import Benchmark
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import numpy as np
from lhstools.utils import discrete_cmap

class Uptake02(Benchmark):
    """Global Uptake 2002-2012"""
    def __init__(self,config_name='config.ini',init=False,*args,**kwargs):
        #Shared parameters as attributes
        Benchmark.__init__(self,config_name,init,*args,**kwargs)
        #FAPAR parameters,
        Benchmark.import_config(self,config_name,'Uptake')
        #Get observations
        self.get_obs()
        #name
        self.name='Uptake02'

    def get_obs(self):
        """set data according to
        IPCC Working Group 1
        """
        self.obs=pd.Series(index=['2002-2011'])

        #From IPCC Ar5
        self.obs['2002-2011']=1.6
        if self.grosscorrect=='True':
            #add correction for Gross Landuse
            self.obs['2002-2011']=self.obs['2002-2011']+0.3471
        self.sigma_obs=1.0



    def calc_stats(self,memberid):
        "Returns stats (error and variance) of a member id"
        #NEW: error of every month
        uptake=self.get_sim(memberid)
        return Benchmark.calc_metric(self,uptake,self.obs,weight=None,sigma_obs=self.sigma_obs)


    def get_sim(self,memberid):
        "Returns the  land emissions from the totc ascii of a member as a pandas dataframe"
        fnmember=self.path2ascii+'trans_'+memberid+'.totc.out'
        totc=Benchmark.read_ascii(self,fnmember)['Total']
        #average over whole simulated period...
        uptake=pd.Series(index=['2002-2011'])
        uptake['2002-2011']=(totc[2011]-totc[2001])/10.
        return uptake

    def plot_hist(self):
        "Pltos a histogram of all the vegetation for all members"
        #Import all the data
        matplotlib.rcParams.update({'font.size': 14})
        uptake=pd.DataFrame(index=self.obs.index)
        for member in self.members:
            uptake[member]=self.get_sim(member)
        # ax=uptake.transpose().hist()
        fig,ax=plt.subplots(figsize=(12,8))
        ax=uptake.transpose().hist(ax=ax, fc='lightblue', histtype='stepfilled', alpha=0.3, normed=True,label='Norm. Histogram',bins=20)
        #Add observations and kde
        uptake.transpose()['2002-2011'].plot(kind='kde',label='KDE',ax=ax[0],bw_method=0.5,color='lightblue')
        ax[0].plot([self.obs['2002-2011'],self.obs['2002-2011']],ax[0].get_ylim(),color='red',linewidth=2,label='IPCC')
        ax[0].axvspan(self.obs['2002-2011']-self.sigma_obs,self.obs['2002-2011']+self.sigma_obs,color='red',alpha=0.3,label=r'IPCC $\pm \sigma$')
        ax[0].legend(loc='upper right',fancybox=True,framealpha=0.8)
        ax[0].set_xlabel('Atmosphere-Land Flux [PgC/yr]')
        return fig,ax
