# ParaFit and PACo analyses
---

Files included in this repo:
- `parafit_script-2024.R`: updated (December 2024) ParaFit script for distance-based host-parasite cophylogenetic analyses
  - loops over multiple (randomized) ParaFit runs
  - adjusts p-values for multiple tests
  - generates simple output of global value, mean p-value, and link p-values
  - plots tanglegram
  - produces a log file to cross check results
  - runs quite a bit faster than the 2017 code
    
- `fake-data-2024.R`: same as above, but modified for fake dataset (data files below)
  - `fake-host.nexus`: fake host dataset
  - `fake-parasite.nexus`: fake parasite dataset
  - `fake-matrix.txt`: fake host-parasite matrix

- `COMING SOON`: PACo analyses
