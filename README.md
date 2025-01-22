# SIMPLICITY 

This is a quick-start guide. For more detailed information please refer to the [wiki](https://gitlab.com/combio/simplicity-public/-/wikis/home).

## Description



## Getting Started
To use SIMPLICITY for running simulations, you will need to setup a python environment with all the required dependencies on your computer and download the code from this repository. In the following sections, we will guide you step by step to install, download and run the code in this repository. We highly reccomend that before running the simulations, you get familiar with the theoretical framework behind this work by referring to our original manuscript.

### Prerequisites
To run SIMPLICITY you need Python3 installed on your computer. After successfully installing Python, you can manually install the dependencies below using pip or Anaconda package manager. Alternatively, you can use the [yml file](https://gitlab.com/combio/simplicity-public/-/blob/main/simplicity_env.yml) provided in the code repository to install all required packages in one step. In Anaconda, you can do so by:

```
conda env create -f simplicity_env.yml
```

For more information on how to use conda environments, please refer to the [official guide](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

conda channels:

- defaults
- conda-forge
- bioconda

dependencies:

- anytree
- biopython
- matplotlib
- numpy
- pandas
- scipy
- scikit-learn
- scikit-posthocs

### Installation

SIMPLICITY do not require any installation. Once you set up the python environment, simply download the code from the repository:

```bash
mkdir path/to/my/folder
git clone https://gitlab.com/combio/simplicity-public
```

## Usage

Before being able to run the code, you need to activate the environment you set up in the previous steps:

```
conda activate simplicity
```

Also make sure to move into the correct folder before running any script or code:

```bash
cd path/to/my/folder/simplicity-public
```

Now you are ready to start using SIMPLICITY! If you want to run simulations from a Python interpreter, you will need to either use one of the provided tuned_parameters json files or to create your own. This file contains all the model parameters you need to specify to run a simulation. Below is an example file you can use as reference when setting up your simulation:

```
{
    "population_size": 1000, # size of the simulated population (number of infected + susceptibles at any time point)
    "infected_individuals_at_start": 10, # number of infected individuals at start
    "R": 1.3, # relative reproduction number of the virus 
    "pT": 11.41, # average infectious time of infected individuals. DO NOT CHANGE unless you know what you are doing.
    "diagnosis rate": 0.0005, # diagnosis rate parameter. Please refer to the corresponding table to see the values that correspond to the desired individuals diagnosis rates.
    "IH virus emergence rate": 0.085, # rate of emergence of intra-host lineages
    "evolutionary rate": 0.0023, # model evolutionary rate. You need to fit this parameter for it to correspond to any desired observed evolutionary rate u
    "t_0": 0, # starting time of the simulation in days
    "final_time": 1095, # final time of the simulation in days (simulation lenght)
    "max_runtime": 100000, # max runtime of simulation in seconds
    "phenotype_model": "distance_from_wt", # phenotype model selection (the other option is immune_waning)
    "sequencing_rate": 0.10, # sequencing rate of diagnosed individuals
    "seed": null, # random seed
    "F": 1.25 # extrande upper bound B factor. DO NOT CHANGE unless you know what you are doing.
}
```

One you have your parameters_file.json ready, you can run a simulation as follows:

```python
import Simplicity
import Population

# create simplicity instance from parameters file
simplicity = Simplicity.SIMPLICITY_simulation(parameters_file)
# create population instance from parameters file
population = Population.create_population(parameters_file)
        
# run SIMPLICITY simulation
simplicity.run(population)
# plot results
simplicity.plot()

# build infection tree    
simplicity.infection_tree()

# build phylogenetic tree    
simplicity.phylogenetic_tree()
```

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.
