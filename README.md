# Chain Files with Isomer Branching ratios

This repo is my attempt at creating an OpenMC depletion chain file with isomer branching ratios for more than (n,gamma) reactions
The scripts have only been tested on the JENDL-5 data library at 14MeV and the process as a whole is still slightly crude.

Q values for isomer producing reactions currently default to the Q value of the ground state producing reaction, I may change this.

ENDF Files are parsed with ENDF-parserpy, data is currently extracted at exactly 14MeV so there is a risk of hitting a resonance peak causing inaccuracy. May update to average over a more reasonable fusion energy spectrum.

Handles branching ratio data stored in both MF=9 and MF=10 subsections, extracts product LFS data from MF=8 then scans decay library to match product ZA number and LFS energy to ZA number and isomeric energy with slight tolerance.
Non matches are sent to ground state.






