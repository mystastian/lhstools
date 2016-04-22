import collections
import xarray as xr
import numpy as np
import pandas as pd
import ast
#Ensure py2&3 compatibiility...
try:
    from configparser import SafeConfigParser
except:
    from ConfigParser import SafeConfigParser

class Benchmark(object):
    """"benchmark parent class"""
    def __init__(self,config_name='config.ini',init=False):
        "Initialization routine, if init=False tries to load certain files from existing files"
        #Import all parameters from the configuration file
        self.import_config(config_name,['General','Advanced'],listnames=['parnames','ignoremembers','lunames'],intnames=['nsamples','setid_start'])
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
        #Create Pandas df to store the parameter values for each member
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
            tempframe=pd.read_csv(self.path2run+member+'.lpj.parameter',delim_whitespace=True,index_col=0,usecols=[0,1],comment='#')
            for paraname in paranames:
                self.para.loc[member,paraname]=float(tempframe.ix[paraname])
        return self.para



    def calc_skill(self):
        "Calculates the skill of all members"
        #error array
        error=pd.Series(index=self.members)
        var=pd.Series(index=self.members)
        for member in self.members:
            error[member],var[member]=self.calc_stats(member)
            if self.verbose:
                print("Member: {} , {}: {}, var:{}".format(member,self.metric,error[member],var[member]))
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
        return pd.read_csv(fn,delim_whitespace=True,header=None,index_col=0,names=names)

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

    def calc_metric(self,model,observation,weight):
        "Dispatch method for error calculation"
        if self.metric=='RMSE':
            return self.calc_rmse(model,observation,weight)
        elif self.metric=='MSE':
            return self.calc_mse(model,observation,weight)
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

    def calc_mse(self,model,observation,weight):
        "Returns the root mean square and variance of model-observation"
        #Case data is provided as xarray
        if isinstance(model,xr.core.dataarray.DataArray) and isinstance(observation,xr.core.dataararay.DataArray):
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
            if weight is None:
                error=model-observation
                se=error**2
                mse=se.mean()
                var=error.var()
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
        pass
