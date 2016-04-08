from ConfigParser import SafeConfigParser
import collections
import xarray as xr

class Benchmark(object):
    """"benchmark parent class"""
    def __init__(self,config_name):
        #Import all parameters from the configuration file
        self.import_config(config_name,['General','Advanced'],intnames=['nsamples','setid_start'])
        self.members=[self.ensembleid+str(i).zfill(4)+self.lhsphase for i in range(self.setid_start,self.nsamples)]


    def import_config(self,config_name,sections,intnames=[]):
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

    def calc_metric(self,model,observation,weight=None):
        "Dispatch method for error calculation"
        if self.metric=='rmse':
            return self.calc_rmse(model,observation,weight)
        else:
            raise ValueError('{} is not a valid metric'.format(self.metric))

    def calc_rmse(self,model,observation,weight):
        "Returns the root mean square of model-observation"
        #Case data is provided as xarray
        if (isinstance(model,xr.core.dataarray.DataArray) and isinstance(observation,xr.core.dataarray.DataArray)):
            return (((model-observation)**2).sum())**0.5
        else:
            raise ValueError('RMSE for {} and {} not implemented'.format(type(model),type(observation)))

    def calc_error(self):
        "Dummy method, should be overloaded by children"
        raise ValueError("method calc_error needs to be overloaded in children")


    def calc_skill(self):
        "Calculates the skill of all members"
        #error array
        error={}
        for member in self.members:
            error[member]=self.calc_error(member)
        return error

