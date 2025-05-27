# SIMPLICITY

SIMPLICITY stands for Stochastic sIMulation of SARS-CoV-2 sPreading and evoLutIon aCcountIng for wiThin-host dYnamics. It is the implementation of a multi-layer, agent based, stochastic mathematical model to study viral spread and evolution and the influence of within host dynamics on population level transmission and evolutionary dynamics. For more detailed info, please refer to the official publication.


[![Documentation Status](https://readthedocs.org/projects/simplicity/badge/?version=latest)](https://simplicity.readthedocs.io/en/latest/)
[![License](https://img.shields.io/github/license/PietroGx/SIMPLICITY)](LICENSE)
[![Python Version](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/)

---

## 📖 Documentation

Full documentation is available at:  
👉 [https://simplicity.readthedocs.io](https://simplicity.readthedocs.io)

---

## Getting Started

### Installation

> 🐍 Requires **Python 3.10 or later**

Clone the repository and install dependencies:

```bash
git clone https://github.com/PietroGx/SIMPLICITY.git
cd SIMPLICITY
pip install -r docs/requirements.txt
```
or use the simplicity.yml file to create the conda environment.



---
## SIMPLICITY Usage and Output

### Usage

Before running a simulation, it's important to understand the following terminology:

| Term               | Definition                                                                 |
|--------------------|----------------------------------------------------------------------------|
| **Experiment**     | A group of simulations with varying parameters used to investigate a specific research or modeling question. |
| **Simulation**     | A set of seeded simulations with shared parameters and different random seeds. |
| **Seeded simulation** | A single simulation with a unique set of parameters and a random seed.   |

To execute an experiment, use the `runme.py` script with the following parameters:

- `experiment_name`: **str**, name of the experiment
- `experiment_settings`: **function**, specifies or loads experiment settings
- `simplicity_runner`: **module**, runner for the experiment
- `plot_trajectory`: **bool**, whether to plot the main simulation output: compartment/fitness/lineage graphs
- `archive_experiment`: **bool**, archive results in a compressed folder

You can also run an example via:

```bash
scripts/experiments/example_experiment_runme.py
```

---

### Model Configuration

Simple changes can be made via `fixture_experiment_settings()`. For more complex changes, edit or generate a parameter `.json` config file using the provided helper script.

---

### Intra-host Model Parameters

The intra-host SARS-CoV-2 dynamics model developed by Van der Toorn
et Al. was designed to evaluate the effectiveness of various testing and
quarantine strategies for managing incoming travelers, contact tracing, and
de-isolation. In our model we use it to specify the individuals
infectiousness and detectability profiles, which depends on the residence
time in a given intra-host model compartment, controlled with the taui
parameters. The user can set their values to obtain different infectiousness
and detectability time windows. In the SIMPLICITY model version
presented here, we only use the parameters from the original publication,
which are fit on german SARS-CoV-2 clinical data. In principle this model
can be adapted, by fitting it to clinical data, to any human pathogenic
virus. It could be interesting, for example, to employ it to model other
respiratory viruses such as influenza. Another possible modification of
intra-host model could be to generate different population types: for
example, if one wanted to model a subpopulation that develops long-covid
or chronic infections, one can modify the intra-host model accordingly.

---

### Intra-host Virus Compartmentalization

In a SIMPLICITY simulation we can simulate the process of virus
compartmentalization, i.e. the establishment of distinct intra-host
dominant lineages that colonize different biological nieces within the host’s
body. The parameter kv (IH virus emergence rate) controls the rate at
which new intra-host dominant lineages are introduced in the individuals.

---

### SIRD Model Parameters


The SIRD model is a classic compartment model consisting of susceptible,
infected, recovered and diagnosed individuals. We assume a fixed
population size (constant number of infected + susceptible individuals)
and well-mixed individuals. The user-set parameters for the model are:

- Population size
- Infected individuals at start
- R
- Diagnosis rate

The infection number R is the simulation parameter that determines
how fast the infection will spread in the population. It controls the average
number of secondary infections caused by each agent in the population.
The diagnosis rate is given in decimal percentage and it is converted
automatically in the corresponding kd value

---

### Evolutionary Model Parameters


The evolutionary model introduces substitutions in the viral genome using
user-defined parameters. Only substitutions are modeled, and the user can
specify the nucleotide substitution rate (NSR), the per-site substitution
rate and a substitution matrix that defines transition and transversion
rates. In our model, we distinguish the NSR from the observed
evolutionary rate (OSR). The NSR denotes the rate at which substitutions
are introduced into the viral genome during the simulation time, while the
OSR is inferred from the sequencing data generated during a simulation
run. Importantly, the OSR does not necessarily match the NSR due to the
effects of population-level transmission and intra-host dynamics. To help
users select an NSR appropriate for their simulation scenario, we provide a
calibration method. During simulations, genomes from diagnosed
individuals are saved at a user-defined sequencing rate. These sequences
are then used to estimate the OSR at the population level via root-to-tip
regression (see Methods). We fit an exponential function to explain the
relationship between NSR and OSR, providing a regression model that
allows users to infer the NSR from the OSR (Figure 3). This allows to
replicate the empirical evolutionary rates observed for SARS-CoV-2 during
the COVID-19 pandemic (in the range of 10−3 − 10−4 nucleotides per
site/per year)

---

### Phenotypic Model Selection

Two models are available:

1. **Distance-from-WT**: Fitness based on average Hamming distance from wild-type
2. **Immune waning**: Accounts for immune memory decay and population exposure

The genotypes of the simulated virus lineages that evolve during a
SIMPLICITY simulation influence the spread dynamics of each lineage.
The phenotypic model assigns a relative fitness score to each viral lineage
based on their genotype and how it differs from what has been circulating
in the population. We implicitly model the host immune system reaction
to viral exposure by assigning fitness score that favors lineages that diverge
from the initial genotype or from the most common circulating lineage over
a time window (implicit immune-escape modeling). We currently have
implemented two versions of the phenotypic model. The first, simpler
version has no parameters: the relative fitness score of an infected
individual is the average of the Hamming distance between their intra-host
lineages and the wild-type (the sequence at the beginning of the
simulation). The second model introduces immune-waning dynamics in the
simulation and its parameters are related to the pharmacokinetics of
SARS-CoV-2 antibodies in human plasma (see Methods). Before running a
SIMPLICITY simulation the user must select which phenotype model to
use.

---

### Model Output

Each seeded simulation generates:

| File Name                        | Description |
|----------------------------------|-------------|
| `final_time.csv`                 | Simulation end time |
| `fitness_trajectory.csv`        | Avg. fitness over time |
| `individuals_data.csv`          | Metadata for each individual |
| `lineage_frequency.csv`         | Frequency of viral lineages |
| `phylogenetic_data.csv`         | Information about lineages |
| `sequencing_data_regression.csv`| Data for estimating OSR |
| `sequencing_data.fasta`         | Simulated sequences |
| `simulation_trajectory.csv`     | SIDR compartment data |

Each SIMPLICITY seeded simulation generates a set of output files
containing the intra-host, population and evolutionary dynamics data. 
These include the compartment trajectories of the SIDR model
(simulation trajectory.csv), individual-level metadata for all infected,
diagnosed, and recovered individuals (individuals data.csv), as well as the
average fitness over time (fitness trajectory.csv) and lineage frequencies
(lineage frequency.csv). Sequencing data — sampled at user-defined rates
from diagnosed individuals — is stored in FASTA format (sequencing
data.fasta), and as regression-formatted data (sequencing data
regression.csv), used to estimate the observed substitution rate. Our model
also provides ground truth transmission and phylogenetic data that can be
used to reconstruct the corresponding tree. 

### Full Model Parameters Table

| Parameter Name              | Description |
|----------------------------|-------------|
| `population_size`          | Total population (infected + susceptible) |
| `infected_individuals_at_start` | Number infected at start |
| `final_time`               | Simulation duration in days |
| `tau_3`                    | Time in infectious/symptomatic state |
| `R`                        | Reproductive number |
| `diagnosis_rate`           | % of population diagnosed |
| `IH_virus_emergence_rate`  | `k_v` — emergence rate of intra-host variants |
| `nucleotide_substitution_rate` | Substitution rate per site/year |
| `phenotype_model`          | Model type (`distance from wt`, `immune waning`) |
| `sequencing_rate`          | % of diagnosed individuals sequenced |
| `max_runtime`              | Max time per seeded simulation |
| `seed`                     | Random seed |

## Project Structure

```plaintext
SIMPLICITY/
├── scripts/                  # scripts, experiments and plots (used in official paper)
├── simplicity/               # Core framework modules
│   ├── evolution/            # Evolutionary algorithm logic
│   ├── phenotype/            # Phenotype model (infection fitness)
│   ├── runners/              # Runners (local/slurm)
│   ├── tree/                 # Infection and phylogenetic tree builders
│   └── tuning/               # Parameter tuning 
├── docs/                     # Documentation (Sphinx)
└── tests/                    # Tests
```

---


## Contributing

Contributions are welcome! If you’d like to report a bug or propose a feature:

1. Open an issue
2. Fork the repo
3. Submit a pull request


---

## 📄 License

This project is licensed under the GPL3 licence. See the [LICENSE](LICENSE) file for details.

---

## 👥 Authors

- Pietro Gerletti
- Jean-Baptiste Escudié

---

## 🌟 Acknowledgements

