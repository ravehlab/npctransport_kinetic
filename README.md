# npctransport_kinetic
Synopsis:
Kinetic model of transport through the NPC, used for comparison with e.g. FLIP experiments

Authors:
Kessem Clein;
Barak Raveh barak.raveh@mail.huji.ac.il

Files:
transport_simulation.py - class for simulating transport
slurm/ - files for running on HPC environment using Slurm (tested on Phoenix cluster, HUJI)
MW_stats_by_force.py - script for testing statistics for different molecular weights with and without force

Unit testing:
Run python test_transport_simulation.py to test that the simulation behaves as expected

Last changed after:
Dec 15, 2020
