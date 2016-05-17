import lhstools as lh
import pandas as pd
import matplotlib.pyplot as plt
from IPython import get_ipython
import xarray as xr
import cartopy.crs as ccrs
import numpy as np

#Settings
_skill=True
_plot=False

eid='S08'
cfg='config_'+eid+'.ini'
saveplts='plots/'+eid+'/'

ipython = get_ipython()
ipython.magic('load_ext autoreload')
ipython.magic('autoreload 2')
plt.ion()

f=lh.FAPAR(cfg)
e=lh.EMDI(cfg)
u=lh.Uptake(cfg)
gv=lh.GlobalVeg(cfg)
gs=lh.GlobalSoil(cfg)
flux=lh.FLUXNET(cfg)

Benchmarks=[u,e,gv,gs,flux]

for bench in Benchmarks:
    try:
        bench.load_stats()
        print(bench.name+' stats loaded')
    except:
        print(bench.name+' calculating stats')
        bench.calc_skill()
        bench.save_stats()


if _skill:
    MSE=pd.DataFrame(index=f.stats.index)
    for bench in Benchmarks:
        if bench!=u:
            MSE[bench.name]=bench.stats[bench.name+'.MSE_rel']
    tskill=lh.Benchmark(cfg)
    tskill.stats[tskill.name+'.skill']=np.exp(-0.5*MSE.sum(axis=1))



if _plot:
#Shared plots
#############
    for bench in Benchmarks:
        fig,ax=bench.plot_skill_vs_para()
        fig.savefig(saveplts+bench.name+'_skill_para.png')
        fig,ax=bench.plot_skill()
        fig.savefig(saveplts+bench.name+'_skill.png')
#Specific plots:
################
    for bench in [u,gv,gs]:
        fig,ax=bench.plot_hist()
        fig.savefig(saveplts+bench.name+'hist.png')
#Bonus Plots

#Soilcarbon
    if gs in Benchmarks:
        soil=pd.DataFrame(index=gs.obs.index)
        for member in gs.members:
            soil[member]=gs.get_sim(member)
        fig,ax=plt.subplots(2,1)
        ax[0].scatter(gs.para['ksoil_tune'],soil)
        ax[0].set_xlabel('ksoil_tune')
        ax[0].set_ylabel('soil carbon')
        try:
            ax[1].scatter(gs.para['E0_hr'],soil)
            ax[1].set_xlabel('E0_hr')
            ax[1].set_ylabel('soil carbon')
        except:
            pass
        fig.savefig(saveplts+gs.name+'_para_soil.png')




