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
        """get data (Land uptake) according to
        Contributions to accelerating atmospheric CO2 growth
        from economic activity, carbon intensity, and efficiency of natural sinks
        Canadell et. al. 2007
        """
        #LUC-Land uptake in PgC/yr
        self.obs=pd.Series(index=['59-10','70-99','90-99','00-06'])
        #LUC + Land uptake
        self.obs['59-10']=-1.5+1.9
        self.obs['70-99']=-1.5+2.0
        self.obs['90-99']=-1.6+2.7
        self.obs['00-06']=-1.5+2.8


    def calc_stats(self,memberid):
        "Returns stats (error and variance) of a member id"
        #NEW: error of every month
        uptake=self.get_sim(memberid)
        return Benchmark.calc_metric(self,uptake,self.obs,weight=None)


    def get_sim(self,memberid):
        "Returns the  land emissions from the totc ascii of a member as a pandas dataframe"
        fnmember=self.path2ascii+'trans_'+memberid+'.totc.out'
        totc=Benchmark.read_ascii(self,fnmember)['Total']
        uptake=pd.Series(index=['59-10','70-99','90-99','00-06'])
        uptake['59-10']=(totc[2010]-totc[1959])/(2010-1959+1)
        uptake['70-99']=(totc[1999]-totc[1970])/(1999-1979+1)
        uptake['90-99']=(totc[1999]-totc[1990])/(1999-1990+1)
        uptake['00-06']=(totc[2006]-totc[2000])/(2006-2000+1)
        return uptake




