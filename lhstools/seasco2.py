from .benchmark import Benchmark
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import numpy as np
from subprocess import call


class SeasCO2(Benchmark):
    """Seasonal CO2 Benchmark using the Globalview-CO2 product and TM2 for atmospheric transport."""
    def __init__(self,config_name='config.ini',init=False,init_tm2=True,*args,**kwargs):
        #Shared parameters as attributes
        super(SeasCO2,self).__init__(config_name,init,*args,**kwargs)
        #SeasCO2 parameters,
        super(SeasCO2,self).import_config(config_name,'SeasCO2')
        #Run TM2 only if init and init_tm2
        if init==False:
            init_tm2=False
        if init_tm2:
            self.run_tm2()
        #Open Ctrl Simulation (After Tracermodel was applied to it [Automatic: ToDo]) for station names
        ds=xr.open_dataset(self.fnctrl[:-4]+'.co2_tm2.nc')
        self.stations=ds.station_name.values
        #Days in Month used for weighted averaging
        self.monthwgt=[31,28,31,30,31,30,31,31,30,31,30,31]
        #Define Timeaxis
        self.timeaxis=pd.date_range(self.startyear+'-01-01',self.endyear+'-01-01',freq='m')

        #Get Observations:
        self.obsmean=pd.DataFrame(index=range(1,13))
        self.obsstd=self.obsmean.copy()
        for station in self.stations:
            obs=pd.read_csv(self.path2obs+station+'_seas.co2',delim_whitespace=True,skiprows=19,header=None,index_col=0)
            self.obsmean[station]=obs[1]
            self.obsstd[station]=obs[2]

        self.name='SeasCO2'

    def run_tm2(self):
        "Executes the ncl script for every member, to calculate the atm. conc. at the 9 stations"
        for memberid in self.members:
            print("Calculation TM2 for member: "+memberid)
            cmd='ncl {tm2path} \'runname="{runname}"\' \'path="{path2cdf}"\' startyear={startyear} endyear={endyear}'.format(
                runname=memberid,path2cdf=self.path2cdf,startyear=self.startyear,endyear=self.endyear,tm2path=self.path2tm2)
            call(cmd,shell=True)

    def calc_stats(self,memberid):
        "Returns stats (error and variance) of a member id"
        #NEW: error of every month
        simmean,_,_,_,_,_=self.get_sim(memberid)
        return super(SeasCO2,self).calc_metric(simmean.stack(),self.obsmean.stack(),None)


    def get_sim(self,memberid):
        "returns tm2 output"
        fnmember=self.path2cdf+memberid+'.co2_tm2.nc'
        #Get Dataset
        ds=xr.open_dataset(fnmember)
        ds.STATION.values=self.stations
        ds.TIME.values=self.timeaxis
        #Setup output
        simmean=pd.DataFrame(index=range(1,13))
        simmean.index.name='month'
        simstd=simmean.copy()
        simstd_land=simmean.copy()
        simstd_ocean=simmean.copy()
        simmean_land=simmean.copy()
        simmean_ocean=simmean.copy()
        for station in self.stations:
            series=ds.sel(STATION=station).conc.to_series()
            series_land=ds.sel(STATION=station).conc_land.to_series()
            series_ocean=ds.sel(STATION=station).conc_ocn.to_series()
            #Set annual offset with days in month weight to zero
            aoffset=series.groupby(series.index.year).apply(lambda x: np.average(x,weights=self.monthwgt))
            aoffset_land=series_land.groupby(series.index.year).apply(lambda x: np.average(x,weights=self.monthwgt))
            aoffset_ocean=series_ocean.groupby(series.index.year).apply(lambda x: np.average(x,weights=self.monthwgt))
            corrseries=series.copy()
            corrseries_land=series.copy()
            corrseries_ocean=series.copy()
            #Correct Series with annual offset
            for y,year in enumerate(range(1990,2010)):
                corrseries.iloc[(12*y):(12*y+12)]=series.iloc[(12*y):(12*y+12)]-aoffset.iloc[y]
                corrseries_land.iloc[(12*y):(12*y+12)]=series_land.iloc[(12*y):(12*y+12)]-aoffset_land.iloc[y]
                corrseries_ocean.iloc[(12*y):(12*y+12)]=series_ocean.iloc[(12*y):(12*y+12)]-aoffset_ocean.iloc[y]

            #Get Seasonality
            simmean[station]=corrseries.groupby(corrseries.index.month).mean()
            simstd[station]=corrseries.groupby(corrseries.index.month).std()
            simmean_land[station]=corrseries_land.groupby(corrseries.index.month).mean()
            simstd_land[station]=corrseries_land.groupby(corrseries.index.month).std()
            simmean_ocean[station]=corrseries_ocean.groupby(corrseries.index.month).mean()
            simstd_ocean[station]=corrseries_ocean.groupby(corrseries.index.month).std()

        return simmean,simstd,simmean_land,simstd_land,simmean_ocean,simstd_ocean


    def plot_member(self,memberid):
        "Plots given member and observed CO2"
        simmean,simstd,simmean_land,simstd_land,simmean_ocean,simstd_ocean=self.get_sim(memberid)
        fig,axs=plt.subplots(3,3,figsize=(16,16))
        axs=axs.flatten()
        for i,station in enumerate(self.stations):
            self.obsmean[station].plot(fmt='o',yerr=self.obsstd[station],color='black',ax=axs[i])
            axs[i].set_xlim(0.5,12.5)
            axs[i].set_ylim(-12.,9.)
            axs[i].set_title('Station: '+station)

            simmean[station].plot(ax=axs[i],color='r')
            axs[i].fill_between(simmean[station].index, simmean[station]-simstd[station], simmean[station]+simstd[station], color='r', alpha=0.2)
#land only
            simmean_land[station].plot(ax=axs[i],color='g')
            axs[i].fill_between(simmean_land[station].index, simmean_land[station]-simstd_land[station], simmean_land[station]+simstd_land[station], color='g', alpha=0.2)
#ocean only
            simmean_ocean[station].plot(ax=axs[i],color='b')
            axs[i].fill_between(simmean_ocean[station].index, simmean_ocean[station]-simstd_ocean[station], simmean_ocean[station]+simstd_ocean[station], color='b', alpha=0.2)
        plt.tight_layout()
        return fig



