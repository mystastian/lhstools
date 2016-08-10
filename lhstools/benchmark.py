import collections
import xarray as xr
import numpy as np
import pandas as pd
import ast
import os.path
import matplotlib.pyplot as plt
#Ensure py2&3 compatibiility...
try:
    from configparser import SafeConfigParser
except:
    from ConfigParser import SafeConfigParser

class Benchmark(object):
    """"benchmark parent class"""
    def __init__(self,config_name='config_S12.ini',init=False):
        "Initialization routine, if init=False tries to load certain files from existing files"
        #Import all parameters from the configuration file
        self.import_config(config_name,['General','Advanced'],listnames=['parnames','parlog','ignoremembers','lunames'],intnames=['nsamples','setid_start'])
        #Get all member names and ignore the ones specified in $ignoremembers
        self.members=[self.ensembleid+str(i).zfill(4)+self.lhsphase for i in range(self.setid_start,self.nsamples)]
        for imember in self.ignoremembers:
            try:
                self.members.remove(self.ensembleid+str(imember).zfill(4)+self.lhsphase)
            except:
                print('WARNING: Failed to remove member {} from self.members'.format(str(imember)))
        self.name='Benchmark'
        #Create Pandas dataframe to store statistics about ensemble
        self.stats=pd.DataFrame(index=self.members)
        self.stats.index.name='member'
        #Check if saved parameter file exists, if not override init
        if not init:
            if not os.path.isfile(self.path2output+self.ensembleid+'_parameters.csv'):
                init=True
        if init:
            self.para=pd.DataFrame(index=self.members)
            self.para.index.name='member'
            self.import_para(self.parnames)
            self.save_para()

        else:
            self.load_para()


    def import_config(self,config_name,sections,listnames=[],intnames=[]):
        "Imports the parameters of a section specified in the configuration file as class attributes"
        #If a single section/name is provided convert to list
        if isinstance(sections,str):
            sections=[sections]
        if isinstance(intnames,str):
            intnames=[intnames]
        parser = SafeConfigParser()
        parser.optionxform = str  # make option names case sensitive
        found = parser.read(config_name)
        if not found:
            raise ValueError('No config file found!')
        for name in sections:
            if not name in parser.sections():
                raise ValueError('Section {sec} not found in {config}'.format(sec=name,config=config_name))
            self.__dict__.update(parser.items(name))
        #Convert strings to integer where appropriate
        for name in intnames:
            if not (hasattr(self,name)):
                raise ValueError('Parameter {par} not found in {config}'.format(par=name,config=config_name))
            setattr(self,name,int(getattr(self,name)))
        for name in listnames:
            if not (hasattr(self,name)):
                raise ValueError('Parameter {par} not found in {config}'.format(par=name,config=config_name))
            setattr(self,name,ast.literal_eval(getattr(self,name)))

    def import_para(self,paranames):
        "Imports the parameter values from the parameter files of the ensemble"
        for member in self.members:
            print member
            tempframe=pd.read_csv(self.path2run+member+'.lpj.parameter',delim_whitespace=True,index_col=0,usecols=[0,1],comment='#')
            for paraname in paranames:
                self.para.loc[member,paraname]=float(tempframe.ix[paraname])
        return self.para



    def calc_skill(self):
        "Calculates the skill of all members"
        #error array
        print len(self.members)
        error=pd.Series(index=self.members)
        var=pd.Series(index=self.members)
        for member in self.members:
            error[member],var[member]=self.calc_stats(member)
            if self.verbose:
                print("{} Member: {} , {}: {}, var:{}".format(self.name,member,self.metric,error[member],var[member]))
        #Add error to dataframe
        self.stats[self.name+'.'+self.metric]=error
        self.stats[self.name+'.var']=var
        bestid=error.argmin()
        rel_error=error/var[bestid]
        self.stats[self.name+'.'+self.metric+'_rel']=rel_error
        self.stats[self.name+'.skill']=np.exp(-0.5*rel_error)
        #Rank the members
        self.stats[self.name+'.rank']=self.stats[self.name+'.skill'].rank(ascending=False)
        return self.stats

#Saving and loading

    def read_ascii(self,fn):
        "Reads a transient ascii and returns it as a Pandas Dataframe"
        names=['Total']+self.lunames
        out=pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0,names=names)
        #Offset 1 year to account for "wrong" lpx output (1800.5 is rounded to 1801)
        out.index=out.index-1
        return out

    def save_stats(self):
        "saves stats DataFrame to self.ensembleid_self.name.csv"
        self.stats.to_csv(self.path2output+self.ensembleid+'_'+self.name+'.csv')

    def load_stats(self):
        "load stats Dataframe"
        self.stats=pd.read_csv(self.path2output+self.ensembleid+'_'+self.name+'.csv',index_col=0)

    def save_para(self):
        "saves parameter DataFrame to self.ensembleid_parameters.csv"
        self.para.to_csv(self.path2output+self.ensembleid+'_parameters.csv')

    def load_para(self):
        "load para Dataframe"
        self.para=pd.read_csv(self.path2output+self.ensembleid+'_parameters.csv',index_col=0)


#Error Metrics and Variance functions

    def calc_metric(self,model,observation,weight,*args,**kwargs):
        "Dispatch method for error calculation"
        if self.metric=='RMSE':
            return self.calc_rmse(model,observation,weight,*args,**kwargs)
        elif self.metric=='MSE':
            return self.calc_mse(model,observation,weight,*args,**kwargs)
        else:
            raise ValueError('{} is not a valid metric'.format(self.metric))

    def calc_rmse(self,model,observation,weight):
        "Returns the root mean square and variance of model-observation"
        #Case data is provided as xarray
        if isinstance(model,xr.core.dataarray.DataArray):
            se=(observation-model)**2*weight
            mse=np.float(se.sum()/weight.sum())
            # rmse=np.float(np.sqrt(se.sum()/weight.sum()))
            print("Variance is currently not weighted in RMSE")
            var=np.float((model-observation).var())
            # Return rmse and sqrt(var)
            return np.sqrt(mse),np.sqrt(var)
        else:
            raise ValueError('RMSE for {} not implemented'.format(type(model)))

    def calc_mse(self,model,observation,weight,sigma_obs=None):
        "Returns the root mean square and variance of model-observation"
        #Case data is provided as xarray
        if isinstance(model,xr.core.dataarray.DataArray) and isinstance(observation,xr.core.dataarray.DataArray):
            if weight is None:
                raise ValueError('Non-Weighted MSE for xarray not implemented yet!')
            else:
                ws=weight.sum()
                error=observation-model
                se=error**2*weight
                mse=np.float(se.sum()/ws)
                meanerror=(error*weight).sum()/ws
                var=np.float( ((error - meanerror)**2 *weight).sum()/ws)
        elif isinstance(model,pd.core.series.Series) and isinstance(observation,pd.core.series.Series):
            if sigma_obs is None:
                if weight is None:
                    error=model-observation
                    se=error**2
                    mse=se.mean()
                    var=error.var()
                else:
                    raise ValueError('Weighted MSE for pandas series not implemented yet!')
            else:
                if weight is None:
                    error=model-observation
                    se=error**2
                    mse=se.mean()
                    var=sigma_obs**2
                else:
                    raise ValueError('Weighted MSE for pandas series not implemented yet!')

        else:
            raise ValueError('MSE for {} not implemented'.format(type(model)))
        return mse,var


    def calc_stats(self,*args):
        "Dummy method, should be overloaded by children"
        raise ValueError("method calc_stats needs to be overloaded in children")

#Plotting

    def plot_skill_vs_para(self):
        "Plots skill of a given target versus all the parameters"
        npars=len(self.parnames)
        fig,axs=plt.subplots((npars-1)/4+1,4,figsize=(16,16))
        axs=axs.flatten()
        plt.suptitle('Ensemble: '+self.ensembleid+', Target:'+self.name+' - Parameters vs Skill',fontsize=20)
        #Delte the axes not used
        if npars%4 != 0:
            for i in range((4-npars)%4,0):
                axs[i].axis('off')

        for i,pname in enumerate(self.parnames):
            if pname in self.parlog:
                axs[i].set_xscale('log')
            axs[i].scatter(self.para[pname],self.stats[self.name+'.skill'],alpha=0.2)
            axs[i].set_xlabel(pname)
            axs[i].set_ylabel('skill [1]')
        #Fix Suptitle
        plt.tight_layout()
        plt.subplots_adjust(top=0.95)
        return fig,axs


    def plot_skill(self):
        "Histogram +KDE of the skill"
        fig,ax=plt.subplots(figsize=(12,8))
        self.stats[self.name+'.skill'].plot(kind='kde',label='No restriction',ax=ax,bw_method=0.5)
        ax.hist(self.stats[self.name+'.skill'].dropna(), fc='lightblue', histtype='stepfilled', alpha=0.3, normed=True)
        ax.set_title('Ensemble: '+self.ensembleid+', Benchmark: '+self.name+' - skill distribution')
        ax.set_xlabel('skill')
        ax.set_xlim([-0.05,1.05])
        return fig,ax

    def plot_best(self):
        "Plot member with highest skill"
        best=self.stats[self.name+'.rank'].argmin()
        try:
            return self.plot_member(best)
        except:
            raise ValueError('plot_member method not available')

    def plot_worst(self):
        "Plot member with lowest skill"
        worst=self.stats[self.name+'.rank'].argmax()
        try:
            return self.plot_member(worst)
        except:
            raise ValueError('plot_member method not available')
