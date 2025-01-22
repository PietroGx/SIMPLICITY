# SIMPLICITY 

This is a quick-start guide. For more detailed information please refer to the [wiki].

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
git clone https://github.com/PietroGx/SIMPLICITY
```

## Usage

Before being able to run the code, you need to activate the environment you set up in the previous steps:

```
conda activate simplicity
```

Also make sure to move into the correct folder before running any script or code:

```bash
cd path/to/my/folder/simplicity
```

Now you are ready to start using SIMPLICITY! If you want to run simulations from a Python interpreter, you will need to either use one of the provided runme.py files or to create your own. Below is an example file you can use as reference when setting up your simulation:

```python
import simplicity.config as config
import simplicity.settings_manager as sm
import simplicity.output_manager as om
import simplicity.runners.serial 
from   simplicity.runners.unit_run import run_seeded_simulation

experiment_name = 'EXAMPLE_RUNME'

## experiment settings (sm.write_settings arguments)
parameters      = {'evolutionary rate': [0.001]}
n_seeds         = 3

def run_experiment(experiment_name, 
                   simplicity_runner,
                   plot_trajectory,
                   archive_experiment=False):
    print('')
    print('##########################################')
    print('')
    # setup experiment files directories
    config.create_directories(experiment_name)
    # set parameters 
    sm.write_settings(parameters, n_seeds)
    # Write experiment settings file
    sm.write_experiment_settings(experiment_name)
    # write simulation parameters files
    sm.read_settings_and_write_simulation_parameters(experiment_name)
    # write seeded simulation parameters files
    sm.write_seeded_simulation_parameters(experiment_name)
    
    # let one of simplicity.runners run each seeded simulation
    simplicity_runner.run_seeded_simulations(experiment_name, 
                                             run_seeded_simulation,
                                             plot_trajectory)
    
    if archive_experiment: 
        om.archive_experiment(experiment_name)
    
    print('')
    print(f'EXPERIMENT {experiment_name} EXECUTED SUCCESSFULLY.')
    print('')
    print('##########################################')

if __name__ == "__main__":
    
    run_experiment(experiment_name,                      
                   simplicity_runner  = simplicity.runners.serial,
                   plot_trajectory = True,
                   archive_experiment = False)
```
The parameters of a SIMPLICITY simulation can be defined by the user in the runme.py file (when specifying the parameters dictionary). Here is an overview of all the parameters and their standard values:
```
STANDARD_VALUES for SIMPLICITY simulation: 
    
    "population_size": 1000               # size of the simulated population (number of infected + susceptibles at any time point)
    "infected_individuals_at_start": 100  # number of infected individuals at start
    "R": 1.5                              # relative reproduction number of the virus 
    "k_d": 0.0055                         # diagnosis rate parameter. Please refer to the corresponding table to see the values that correspond to the desired individuals diagnosis rates.
    "k_v": 0.0085                         # rate of emergence of intra-host lineages
    "e": 0.0017 (evolutionary rate)       # model evolutionary rate. You need to fit this parameter for it to correspond to any desired observed evolutionary rate u
    "final_time": 365*3                   # final time of the simulation in days (simulation lenght)
    "max_runtime": 300                    # max runtime of simulation in seconds
    "phenotype_model": 'immune waning' or 'distance from wt'
    "sequencing_rate": 0.05               # sequencing rate of diagnosed individuals
    "seed": None                          # random seed
    "F": 1.25                             # extrande upper bound B factor. DO NOT CHANGE unless you know what you are doing.

```

If you want to change any, for each parameter, specify a list of values that you would like to use for the simulation. If you want to change more than one parameter at the time, consider that you need to enter the same number of values for each parameter, e.g. :

    par_1 = [value1, value2]
    par_2 = [value3, value4]

This will run a simulation with par_1 = value1 and par_2 = value3, and a simulation with par_1 = value2 and par_2 = value4. 

Each simulation will be repeated n_seeds time with a different random seed.

The set of all simulations performed in a runme.py is what we call an experiment.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.
