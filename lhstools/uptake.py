from .benchmark import Benchmark
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import numpy as np
from lhstools.utils import discrete_cmap

class Uptake(Benchmark):
    """Land uptake (Uptake+LUC emissions benchmark"""
    def __init__(self,config_name='config.ini',init=False,*args,**kwargs):
        #Shared parameters as attributes
        Benchmark.__init__(self,config_name,init,*args,**kwargs)
        #FAPAR parameters,
        Benchmark.import_config(self,config_name,'Uptake')
        #Get observations
        self.get_obs()
        #name
        self.name='uptake'

    def get_obs(self):
        """get data (Land uptake) according to specified source
        """
        #LUC-Land uptake in PgC/yr
        if self.source=='Canadell':
            # Contributions to accelerating atmospheric CO2 growth
            # from economic activity, carbon intensity, and efficiency of natural sinks
            # Canadell et. al. 2007
            self.obs=pd.Series(index=['59-10','70-99','90-99','00-06'])
            #LUC + Land uptake
            self.obs['59-10']=-1.5+1.9
            self.obs['70-99']=-1.5+2.0
            self.obs['90-99']=-1.6+2.7
            self.obs['00-06']=-1.5+2.8
        elif self.source=='IPCC':
            #IPCC 5th Assesment report chapter 6. page 486
            #LUC + Land uptake
            self.obs=pd.Series(index=['1980-1989','1990-1999','2000-2009','2002-2011'])
            self.obs['1980-1989']=0.1
            self.obs['1990-1999']=1.1
            self.obs['2000-2009']=1.5
            self.obs['2002-2011']=1.6
            self.obs['1750-2011']=-30
        else:
            raise ValueError('{} is not a valid uptake data source'.format(self.source))


    def calc_stats(self,memberid):
        "Returns stats (error and variance) of a member id"
        #NEW: error of every month
        uptake=self.get_sim(memberid)
        return Benchmark.calc_metric(self,uptake,self.obs,weight=None)


    def get_sim(self,memberid):
        "Returns the  land emissions from the totc ascii of a member as a pandas dataframe"
        fnmember=self.path2ascii+'trans_'+memberid+'.totc.out'
        totc=Benchmark.read_ascii(self,fnmember)['Total']
        if self.source=='Canadell':
            uptake=pd.Series(index=['59-10','70-99','90-99','00-06'])
            uptake['59-10']=(totc[2010]-totc[1959])/(2010-1959+1)
            uptake['70-99']=(totc[1999]-totc[1970])/(1999-1979+1)
            uptake['90-99']=(totc[1999]-totc[1990])/(1999-1990+1)
            uptake['00-06']=(totc[2006]-totc[2000])/(2006-2000+1)
        elif self.source=='IPCC':
            uptake=pd.Series(index=['1980-1989','1990-1999','2000-2009','2002-2011'])
            uptake['1980-1989']=(totc[1989]-totc[1979])/10.
            uptake['1990-1999']=(totc[1999]-totc[1989])/10.
            uptake['2000-2009']=(totc[2009]-totc[1999])/10.
            uptake['2002-2011']=(totc[2011]-totc[2001])/10.
            uptake['1750-2011']=(totc[2011]-totc[1801])

        return uptake

    def plot_hist(self):
        "Pltos a histogram of all the uptake rates for all members"
        #Import all the data
        uptake=pd.DataFrame(index=self.obs.index)
        for member in self.members:
            uptake[member]=self.get_sim(member)
        # ax=uptake.transpose().hist()
        fig,ax=plt.subplots(figsize=(12,8))
        ax=uptake.transpose().hist(ax=ax, fc='lightblue', histtype='stepfilled', alpha=0.3, normed=True,label='Norm. Histogram')
        plt.suptitle('Ensemble: '+self.ensembleid+' - Land Uptake for 4 Periods')
        #Add observations and kde
        if self.source=='IPCC':
            uptake.transpose()['1980-1989'].plot(kind='kde',label='KDE',ax=ax[0,0],bw_method=0.5)
            ax[0,0].plot([self.obs['1980-1989'],self.obs['1980-1989']],ax[0,0].get_ylim(),color='red',linewidth=3,label='IPCC')
            uptake.transpose()['1990-1999'].plot(kind='kde',label='KDE',ax=ax[1,0],bw_method=0.5)
            ax[1,0].plot([self.obs['1990-1999'],self.obs['1990-1999']],ax[0,0].get_ylim(),color='red',linewidth=3,label='IPCC')
            uptake.transpose()['2000-2009'].plot(kind='kde',label='KDE',ax=ax[0,1],bw_method=0.5)
            ax[0,1].plot([self.obs['2000-2009'],self.obs['2000-2009']],ax[0,0].get_ylim(),color='red',linewidth=3,label='IPCC')
            uptake.transpose()['2002-2011'].plot(kind='kde',label='KDE',ax=ax[1,1],bw_method=0.5)
            ax[1,1].plot([self.obs['2002-2011'],self.obs['2002-2011']],ax[0,0].get_ylim(),color='red',linewidth=3,label='IPCC')
        ax[1,1].legend(loc='upper left',fancybox=True,framealpha=0.8)
        plt.tight_layout()
        return fig,ax

