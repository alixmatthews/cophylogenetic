## ParaFit and PACo analyses

Two common distance-based analyses used to measure the statistical congruence between host and symbiont phylogenies. The end metrics yield an overall "cophylogenetic congruence" between hosts and their symbionts as well as the contribution of individual host-symbiont links to the overall congruence. This script can easily be modified for different systems (just need tree files and a matrix of associations) and should be of use to many researchers in this field! Please contact me if you discover any bugs or have suggestions for improvement.

ParaFit originally derived from Legendre et al. 2002 *Systematic Biology* "A statistical test for host-parasite coevolution"

PACo originally derived from Balbuena et al. 2013 *PLoS One* "PACo: a novel Procrustes application to cophylogenetic analysis" (and, for R, Hutchinson et al. 2017 *Methods in Ecology and Evolution* "paco: implementing Procrustean approach to cophylogeny in R")

---

Files included in this repo:
- `parafit_script-2024.R`: updated (December 2024) ParaFit script
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
