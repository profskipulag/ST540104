import os
import re
import json
import copy
import glob
import cfgrib
import datetime
import numpy as np
import arviz as az
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from tqdm import tqdm
from cartopy.feature import NaturalEarthFeature
from scipy.stats import binned_statistic
# hack to get Stan to work
import nest_asyncio
nest_asyncio.apply()
import stan






class Forecast:

    def __init__(self, emulator, observations):

        self.emulator = emulator

        self.observations = observations

    def get_stan_data(self):

        ds_puff = self.emulator.da_puff

        result = self.observations

        concs = []
        
        for name in result:
            
            local_id, species, variable = name.split("#")
        
            if variable == 'value':
                
                concs.append(result[name])
                
                concs[-1].name = local_id
        
        df_obs  = (
                    xr
                    .merge(concs)
                    .to_dataframe()
                    .unstack()
                    .reset_index()
                    .rename(columns={'level_0':'local_id','endtime':'date',0:'conc'})
                )
        
        df_obs = df_obs[
                            (df_obs['date']>= ds_puff['date'].values[0]) & 
                            (df_obs['date']<=ds_puff['date'].values[-1])
                        ]
        
        local_ids  = ds_puff.local_id.values
        
        dates = ds_puff.date.values
        
        source_starts = ds_puff['source.source_start'].values
        
        sss = pd.to_datetime(ds_puff["date"].values[0])
        
        sss = datetime.datetime(year=sss.year,month=sss.month, day=sss.day, hour=sss.hour)
        
        source_starts = [ sss + datetime.timedelta(hours = int(hour)) for hour in source_starts]
        
        # round to nearest hour
        dates = np.array(np.array(dates, dtype='datetime64[h]'), dtype='datetime64[s]')
        
        dates  = [pd.Timestamp(d) for d in dates]
        
        lookup_local_id = dict(zip(local_ids, range(len(local_ids))))
        
        lookup_dates = dict(zip(dates, range(len(dates))))
        
        lookup_source_starts = dict(zip(source_starts, range(len(source_starts))))
        
        df_obs['i'] = df_obs['local_id'].apply(lambda r: lookup_local_id[r])
        
        df_obs['j'] = df_obs['date'].apply(lambda r: lookup_dates[r])
        
        df_obs = df_obs.dropna().reset_index()
        
        print("FILTER FOR NEGATIVE SO2 CONC DISABLED!!")
        #df_obs  = df_obs[df_obs['conc']>0].reset_index()

        self.stan_data = {
            
                # dimensions of 'eumlator' array
                'N_puf': len(ds_puff['source.source_start']),
                'N_stn': len(ds_puff.local_id),
                'N_hrs': len(ds_puff.date),
                'N_hts': len(ds_puff['source.height_above_vent']),
            
                # The 'emulator' array. Remember order! [N_puf, N_stn, N_hrs, N_hts] 
                'puffs': ds_puff.values,
            
                # offset and scale for the emulator values
                'height_min': ds_puff['source.height_above_vent'].values[0],
                'height_delta': ds_puff['source.height_above_vent'].values[1]-ds_puff['source.height_above_vent'].values[0],
            
                # number of observations
                #'N_oc':len(df_synth_obs),
                'N_oc':len(df_obs),
                # observations
                #'obs_conc':df_synth_obs['conc plus noise'].values, #/1e6, # micro g -> g
                'obs_conc':df_obs['conc'].values/1e6, # micro g -> g,
                #'ij_conc': df_synth_obs[['i','j']].values+1, # REMEMBER +1 because Stan counts from 1!!
                'ij_conc':df_obs[['i','j']].values+1,
            
                
                # limits for parameter values
                'height_lower': 100,#ds_puff['source.height_above_vent'].values[1]+10,
                'height_upper': 1000,# ds_puff['source.height_above_vent'].values[-2]-10,
            
                'height_mu_loc': 600.0,
                'height_mu_scale':500.0,
                
                'height_sigma_loc': 500.0, 
                'height_sigma_scale':500.0,
                
                
                
                'flux_lower': 0.0,
                'flux_upper':200.0,
            
                'flux_mu_loc':65.0,
                'flux_mu_scale':10.0,
            
                'flux_sigma_loc':30.0,
                'flux_sigma_scale':10.0,
                
            
                
                'sigma_conc_lower': 0.0,
                'sigma_conc_upper': 1e-4,
            
                'sigma_conc_loc':1e-6,
                'sigma_conc_scale':1e-5,
                

                
                'POSTERIOR':0
            }

    def load_stan_model(self):
        with open("pyinfer/model.stan") as f:
    
            stan_code = f.readlines()
        
        self.stan_code  = "".join(stan_code)
        
    def run_prior(self):

        self.stan_data['LIKELIHOOD'] = 0
        self.prior_model = stan.build(program_code=self.stan_code, data=self.stan_data)
        self.prior_fit = self.prior_model.sample(num_chains=4, num_samples=500, delta=0.99)
        

    def run_posterior(self):

        self.stan_data['LIKELIHOOD'] = 1
        self.posterior_model = stan.build(program_code=self.stan_code, data=self.stan_data)
        self.posterior_fit = self.posterior_model.sample(num_chains=4, num_samples=500, delta=0.99)

    def get_inference_data(self):

        ds_puff = self.emulator.da_puff
        
        source_start = [ pd.to_datetime(ds_puff["date"].values[0]) + datetime.timedelta(hours = int(hour)) for hour in ds_puff["source.source_start"].values]

        
        self.inference_data = az.from_pystan(
            posterior = self.posterior_fit,
            prior = self.prior_fit,
            prior_model = self.prior_model,
            posterior_model = self.posterior_model,
            posterior_predictive = ["obs_hat"], 
            prior_predictive = ["obs_hat"],
            observed_data = ["obs_conc"],
            coords = {
                "local_id":ds_puff["local_id"].values,
                "date":ds_puff["date"].values,
                "source_start": source_start,
            },
            dims = {
                "conc_matrix_hat":["chain","draw","local_id","date"],
                "conc_matrix":["chain","draw","local_id","date"],
                "height0":["chain","draw","source_start"],
                "height":["chain","draw","source_start"],
                "flux0":["chain","draw","source_start"],
                "flux":["chain","draw","source_start"],
                "puffs_interpolated":["chain","draw","local_id","date"],
                
                "heightb":["chain","draw","source_start"],
                "fluxb":["chain","draw","source_start"],
        
            }
        )
        
    def summarize_prior(self):

        return az.summary(self.inference_data,group='prior')

    def summarize_posterior(self):

        return az.summary(self.inference_data,group='posterior')

    def prior_traceplot(self):
        
        az.plot_trace(self.inference_data.prior,var_names=[
            'height',
            'flux',
            'sigma_conc',
            'height_mu',
            'height_sigma',
            'flux_mu',
            'flux_sigma'
        ])
        
        plt.tight_layout()

    def posterior_traceplot(self):
        
        az.plot_trace(self.inference_data.posterior,var_names=[
            'height',
            'flux',
            'sigma_conc',
            'height_mu',
            'height_sigma',
            'flux_mu',
            'flux_sigma'
        ])
        
        plt.tight_layout()
        


        
    
                
                

        







    


