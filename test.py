import datetime
import numpy as np
import arviz as az
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
from pyinfer import Forecast
import xarray as xr
from pyem import Emulator

# load observations from netcdf
observations = xr.open_dataset("dt5402.nc")


# load emulator from netdcf
emulator = Emulator.from_netcdf("dt5404.nc")


# The Forecast class provided by the package pyinfer encapsulates the Bayesian inference process,
# which uses the probabilistic programming language Stan (https://mc-stan.org/)

# first we create a new object using the SO2 concentration data and the emulator
forecast = Forecast(emulator=emulator, observations=observations)


forecast.load_stan_model()

# The Bayesian model is specified in the file model.stan in the pyinfer directory
# The data has to be prepared before being passed to stan:
forecast.get_stan_data()

# we then run the model for the prior ....
forecast.run_prior()

# ... and the posterior ...
forecast.run_posterior()

# .. and convert the output to an inferencedata object so arviz (specialised
# Bayesian model output visualisation library) can work with it
forecast.get_inference_data()


# we inspect the results. First, the prior ...
print(forecast.summarize_prior())

# ... export to netcdf for ingestion by 
forecast.to_netcdf("dt5405.nc")

