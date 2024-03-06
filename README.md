This repository contains code to analyze preprocessed data and run figures from the following paper:

%ADD PAPER REFERENCE

Minimally preprocessed data, as well as analyzed data and figures from the paper are available on DABI at https://dabi.loni.usc.edu/dsi/BM2ZIVWKBFH8.

Requirements:
- Most scripts were run using MATLAB R2019b (but likely would work on more recent versions).
- analysis/echo_analysis.m was on MATLAB R2022b (requires higher version than R2019b due to the use of MATLAB function pyrun.m, introduced in R2021b).
- analysis/echo_analysis.py runs on Python 3.9.
- The following packages are required:
	- Chronux (http://chronux.org/)
	- FieldTrip (https://github.com/fieldtrip/fieldtrip)
	- fooof (https://github.com/fooof-tools/fooof_mat)
	- CircStats (https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics)
	- export_fig (https://github.com/altmany/export_fig)
	- Violinplot (https://github.com/bastibe/Violinplot-Matlab)
	- venn (https://www.mathworks.com/matlabcentral/fileexchange/22282-venn)
