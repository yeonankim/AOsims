# ISETBio_AO_UW
## ISETBio Simulations of Adaptive Optics (AO) Experiments from [Sabesan Lab (U. Washington)](http://depts.washington.edu/sabaolab/)

This repository contains demos, simulations, and work-in-progress to model AO experiments using [ISETBio](https://github.com/isetbio).
- Defeat of typical ISETBio modeling of human LCA (custom-enabled by David Brainard)
- Specifications of UW setup (SPD and pixel scaling of displays)
- Demonstrations of UW stimuli, which are built into ISETBio "scenes" using the same code we used to populate PsychToolbox buffers
  - Chromatic Gratings ([Neitz et al. (2020) JOSAA](https://doi.org/10.1364/JOSAA.382384))
  - Tiny red/green bars ([Coates et al. (2019) ARVO abstract](https://iovs.arvojournals.org/article.aspx?articleid=2741793))
  
## Contents

Stable:
- `make_optics.m:` Creates two optics objects, one with typical LCA (ISETBio default) and one with LCA defeated (custom). Seems harmless to construct both objects, even if only one is used in a given simulation.
- `t_twoline_scene.m:` Tutorial showing how to create a two-line scene. Plots a typical R/G bar pair optical image, both with and without LCA.
- `t_twoline_excitations.m:` Tutorial showing how to generate cone excitations from a two-line stimulus (without LCA), using built-in `visualizeConeResponses`

Work-in-progress:
- `GaborSceneNoLCA_postNoise.m:` Extensive simulations of chromatic gratings at various spatial frequencies and conditions, which are fed into an SVM to determine simulated detection thresholds. Currently being cleaned-up.
