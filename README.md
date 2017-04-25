# NuRef
Neutron Reflection (NuRef) project supported through the Nuclear Engineering Design Collaborative (NEDC) club based out of the University of California, Berkeley. Replicates the geometry and physics for the High Flux Neutron Generator (HNFG) accelerator through the Geant4 toolkit for use as a diagnostic tool.

Code modified from the Geant4 example source and header files available through https://github.com/Geant4/geant4/tree/master/examples/

Contributers: Yaroslav Kaminskiy, Alexander Blank, Haotian Zeng, Joseph Michael Gordon, Zirui Jiang, Donald Nordwick

Oversight and project management: Robert Nnaemeka Nnamani

# Histogramming
The following histograms are produced for each volume of interest:

- All Neutrons Entering the Volume
- Neutrons Scattered In the Volume
- Energy Absorbed in the Volume

Producing a set of these histograms for each logical volume in a scoring volume store will be introduced in the next release. The flux and absorption histogram binning will be introduced in a future release.

# Build Files
The build files and other necessary files to run the simulation executable are in the Release folder. Any ROOT and HepRApp files from the latest simulation output are stored in the the Simulation Output folder.

# How to Run
At Idle, type the following: /control/execute gps.mac

Note that the full hfng geometry is not read in the read_gdml.mac file. Instead a geometry that consists of few HFNG parts were used. This geometry is named "test".