# ST540204 infer.py
A repository for finding the posterior distribution over Eruption Source Parameters (ESPs, here, plume height and $SO_2$ flux), as well as ground $SO_2$ concentration, in light of observations, for a given atmospheric dispersion emulator. This can be thought of as a kind of Source Term Estimation (STE) or Data Assimilation (DA), resulting in probabilistic estimates of ESPs and probabilistic forecasts of air quality that are then used for calculating concentration exceedence probabilities. Bayesian inference is performed using Stan, which allows empirical relationships (e.g. between plume height and ambient atmospheric conditions) to be included in the model (not yet implemented). This repository consists of `pyinfer`, a package that contains some helper classes, including 

 * `Forecast`, a class that wraps the pre and post processing of data for performing Bayesian inference with Stan, and provides some functions for visuaslisation and export serving using the API in ST540105.

The useage of `Forecast` is illustrated in the Jupyter notebook `notebook.ipynb`, and the script `infer.py` is the application of the library for the purposes of DTC4. 

## To do
 * Add visualisations as functions to `Forecast` class
 * Add options for empirical relationshipos to the Stan model
 * Add function for exporting data to netcdf for use by web API backend
 * Add user specified priors
 
## Package structure


    ST540104/
    ├── environment.yaml      - configures conda environment with required packages
    ├── LICENSE               - GPL 3
    ├── notebook.ipynb        - example use of the pyem package
    ├── pyinfer               - python inference package
    │   ├── __init__.py       - initialises the package
    │   ├── model.stan        - specification of the Bayesian model
    │   └── source.py         - source code for various classes
    ├── README.md             - this file
    ├── .gitignore            - files to be ignored by git
    └── infer.py              - script that calls the package for DTC4


## To download the repository
Clone the repository to your machine

    git clone https://github.com/profskipulag/ST540104.git

You will be asked for your username and password. For the password github now requires a token:
- on github, click yur user icon in the top right corner
- settings -> developer settings -> personal access tokens -> Tokens (classic) -> Generate new token -> Generate new token (classic) 
- enter you authentifcation code
- under note give it a name, click "repo" to select al check boxes, then click generate token
- copy result enter it as password

## To run the jupyter notebook
Create a new conda environment from the environment.yaml file:

    conda env create -f environment.yaml

Activate the environment

    conda activate st540104
    
Launch the notebook server

    jupyter notebook
    
Navigate to the st540104 directory and click the file `notebook.ipynb` to launch it.