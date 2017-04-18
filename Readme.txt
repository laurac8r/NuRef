Neutron Reflection (NuRef) project supported through the Nuclear Engineering Design Collaborative (NEDC) club based out of the University of California, Berkeley. Replicates the geometry and physics for the High Flux Neutron Generator (HNFG) accelerator through the Geant4 toolkit for use as a diagnostic tool.

Code modified from the Geant4 example source and header files available through https://github.com/Geant4/geant4/tree/master/examples/

Contributers: Yaroslav Kaminskiy, Alexander Blank, Haotian Zeng, Joseph Michael Gordon, Zirui Jiang, Donald Nordwick

Oversight and project management: Robert Nnaemeka Nnamani

---------------
~Histogramming~
---------------
It contains the implementation of historamming of particles entering a volume, one where all particles are binned including incident one, the other where only scattered particles are binned. 
The Build folder is where it was built. One can build ones own. But one need to go to the Release folder inside the Build folder and copy the necessary files to the Release folder of one's own build. 

------------
~How to run~
------------
At Idle, type the following: /control/execute run1.mac

Note that the full hfng geometry is not read in the read_gdml file. Instead a geometry that consists of few hfng parts were used. This geometry is named "test".