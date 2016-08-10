from .benchmark import Benchmark
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import numpy as np
from lhstools.utils import discrete_cmap

class TotalUptake(Benchmark):
    """Land uptake (Uptake+LUC emissions benchmark over whole industrial period"""
    def __init__(self,config_name='config.ini',init=False,*args,**kwargs):
        #Shared parameters as attributes
        Benchmark.__init__(self,config_name,init,*args,**kwargs)
        #FAPAR parameters,
        Benchmark.import_config(self,config_name,'Uptake')
        #Get observations
        self.get_obs()
        #name
        self.name='totaluptake'

    def get_obs(self):
        """get data (Land uptake) according to specified source
        """
        self.obs=pd.Series(index=['1750-2009'])
        self.obs['1750-2011']=-30


    def calc_stats(self,memberid):
        "Returns stats (error and variance) of a member id"
        #NEW: error of every month
        uptake=self.get_sim(memberid)
        return Benchmark.calc_metric(self,uptake,self.obs,weight=None,sigma_obs=45)


    def get_sim(self,memberid):
        "Returns the  land emissions from the totc ascii of a member as a pandas dataframe"
        fnmember=self.path2ascii+'trans_'+memberid+'.totc.out'
        totc=Benchmark.read_ascii(self,fnmember)['Total']
        uptake=pd.Series(index=['1750-2009'])
        #add the period 1801-1851 twice to get an estimate for 1750-1800
        uptake['1750-2011']=(totc[2011]-totc[1801])+totc[1851]-totc[1801]

        return uptake

    def plot_hist(self):
        "Pltos a histogram of all the uptake rates for all members"
        #Import all the data
        uptake=pd.DataFrame(index=self.obs.index)
        for member in self.members:
            uptake[member]=self.get_sim(member)
        # ax=uptake.transpose().hist()
        fig,ax=plt.subplots(2,2,figsize=(12,8))
        # plt.suptitle('Ensemble: '+self.ensembleid+' - Land Uptake for 4 Periods')
        #Add observations and kde
        if self.source=='IPCC':
            uptake.transpose()['1980-1989'].hist(ax=ax[0,0], fc='lightblue', histtype='stepfilled', alpha=0.3, normed=True,label='Norm. Histogram')
            uptake.transpose()['1980-1989'].plot(kind='kde',label='KDE',ax=ax[0,0],bw_method=0.5)

            ax[0,0].set_title('Uptake 1980-1989')
            ax[0,0].set_xlabel(r'Atm.-Land flux [PgC yr^${-1}$]')
            ax[0,0].plot([self.obs['1980-1989'],self.obs['1980-1989']],ax[0,0].get_ylim(),color='red',linewidth=3,label='IPCC')
            uptake.transpose()['1990-1999'].hist(ax=ax[1,0], fc='lightblue', histtype='stepfilled', alpha=0.3, normed=True,label='Norm. Histogram')
            uptake.transpose()['1990-1999'].plot(kind='kde',label='KDE',ax=ax[1,0],bw_method=0.5)
            ax[1,0].plot([self.obs['1990-1999'],self.obs['1990-1999']],ax[0,0].get_ylim(),color='red',linewidth=3,label='IPCC')
            ax[1,0].set_title('Uptake 1990-1999')
            ax[1,0].set_xlabel(r'Atm.-Land flux [PgC yr^${-1}$]')
            uptake.transpose()['2000-2009'].hist(ax=ax[0,1], fc='lightblue', histtype='stepfilled', alpha=0.3, normed=True,label='Norm. Histogram')
            uptake.transpose()['2000-2009'].plot(kind='kde',label='KDE',ax=ax[0,1],bw_method=0.5)
            ax[0,1].plot([self.obs['2000-2009'],self.obs['2000-2009']],ax[0,0].get_ylim(),color='red',linewidth=3,label='IPCC')
            ax[0,1].set_title('Uptake 2000-2009')
            ax[0,1].set_xlabel(r'Atm.-Land flux [PgC yr^${-1}$]')
            uptake.transpose()['2002-2011'].hist(ax=ax[1,1], fc='lightblue', histtype='stepfilled', alpha=0.3, normed=True,label='Norm. Histogram')
            uptake.transpose()['2002-2011'].plot(kind='kde',label='KDE',ax=ax[1,1],bw_method=0.5)
            ax[1,1].plot([self.obs['2002-2011'],self.obs['2002-2011']],ax[0,0].get_ylim(),color='red',linewidth=3,label='IPCC')
            ax[1,1].set_title('Uptake 2002-2011')
            ax[1,1].set_xlabel(r'Atm.-Land flux [PgC yr^${-1}$]')
        ax[1,1].legend(loc='upper left',fancybox=True,framealpha=0.8)
        plt.tight_layout()
        return fig,ax

    def plot_member(self,memberid):
        fnmember=self.path2ascii+'trans_'+memberid+'.totc.out'
        totc=Benchmark.read_ascii(self,fnmember)['Total']
        fig,ax=plt.subplots(figsize=(12,8))
        (totc-totc.ix[1801]).plot()
        ax.set_title('Land Atmosphere Flux, '+memberid)




