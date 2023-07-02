# Recreate the results of M Giraud et al. (2023)

Bellow we describe the steps to recreate the results of M.Giraud et al (2023) about the developments of CPlantBox 2.0:

## Setup the files

1) Download the files in a linux environment.
2) make sure you have all the requirements installed by running the file "DUMUX/checkRequirements.py"
## Recreate the main results
1) run the three simulations presented in the paper (baseline, ealry dry spell, late dry spell), by running:
```
    cd DUMUX/dumux-rosi/python/coupled
    python3 runSimulation.py baseline
    python3 runSimulation.py ealryDry
    python3 runSimulation.py lateDry
```
2) once the simulations have run, plot the results by running all the .R files in the folder "runSimulation.py".
## Recreate the sensitivty analysis results
1) run the sensitivity analsisys for the xylem and phloem modules:
```
    cd DUMUX/dumux-rosi/python/coupled
    python3 sobol_phloem.py 18 dry
    python3 sobol_phloem.py  11 dry
    python3 sobol_phloem.py 18 wet
    python3 sobol_phloem.py  11 dry
    python3 sobol_xylem.py 18 dry
    python3 sobol_xylem.py  11 dry
    python3 sobol_xylem.py 18 wet
    python3 sobol_xylem.py  11 dry
```
Alternativelly, adapt and use the batch files int he "bashFilesCluster" folder.

2) Plot the results by going through the jupyter notebook "plotResults/plotSobolOutput.ipynb"
