# PIE-Model-Experiments
This repository provides a guide and the necessary resources for reproducing the experimental results of the publication: "A Physics-Informed Eikonal Model for Simulating Arrhythmias in the Human Heart in Real-Time".

## Description
The `setups` directory contains the six experimental setups from the study:
* A1_anatomical: an anatomical reentry is simulated by pacing below an isthmus-shaped scar using an S1-S2 protocol.
* A2_functional: a functional reentry is induced by stimulating a patch of tissue that overlaps with the waveback of a preceding planar wave front.
* A3_wholeheart: a scar-induced ventricular tachycardia is simulated in an anatomically detailed human whole heart model. 
* B1_restitution: the left-hand face of a strip of healthy tissue is paced with a dynamic pacing protocol to investigate APD and CV restitution properties.
* B2_curvature: concave and convex wave fronts are induced in healthy tissue to investigate curvature properties.
* B3_diffusion: an interface between healthy and borderzone tissue types is paced at the bottom to investigate diffusion properties.

Each setup consists of a simulation plan in `.json` format and directories with mesh files in 250um and 1000um resolution stored in both binary carp and VTK formats. The binary carp meshes can be visualized in [NumeriCor Studio](https://numericor.at/rlb/wordpress/products/#ProdStudio) while the VTK meshes can be visualized in [ParaView](https://www.paraview.org/). The human whole-heart experiment was derived from [1] and additionally contains a directory with the lead-field data required for ECG reconstruction. 

Note: The meshes for each setup, except the `A3_wholeheart`, are generated at runtime if not present. The complete dataset including the calibration data, the generated meshes, the anonymized human whole-heart model with lead-field data and rendered videos of each experiment are available via [Zenodo](https://doi.org/10.5281/zenodo.17198150) after the embargo is lifted.

The `calibration` directory contains files that represent captured physiological properties aquired through ForCEPSS [2] which are referenced by the simulation plans. 

The `scripts` directory contains a collection of python scripts to generate, process or visualize the experimental data:
* `meshgen.py`: performs automatic mesh generation for each experiment in case of absence.
* `visualize-calibration.py`: visualizes the calibration results.
* `apply_tstart_offset.py`: corrects the timing offset of LAT files caused by prepacing.
* `error-density.py`: plots KDEs of LAT, CV and LRT error density maps for a given setup.
* `phase_singularitiy_tracking.py`: performs phase singularity tracking and visualization.
* `ecg_comparison.py`: visualizes computed ECGs.

The `eval.sh` bash file automatically runs all operations for reproducing the experimental data.

The `clean.sh` bash file cleans generated data within this repository.

## Prequesites
A [python3](https://www.python.org) environment is required to reproduce the experimental data of the PIE model. We recommend to use [Miniconda3](https://www.anaconda.com/docs/getting-started/miniconda/install) to create a suitable environment on your machine: 

```bash
conda create --name pie-env --file requirements.txt
conda activate pie-env
```

The reference data is generated using the [openCARP](https://opencarp.org/) simulation framework which also includes [mesher](https://git.opencarp.org/openCARP/openCARP/-/tree/master/tools/mesher) and [meshtool](https://bitbucket.org/aneic/meshtool/src/master/) required for mesh generation. Documentation and instructions for installing openCARP can be found [here](https://opencarp.org/download/installation). Alternatively, all neccessary prequesites for reproducing the experimental data are available within the [openCARP Docker Image](https://opencarp.org/download/installation#installation-of-opencarp-).

## Setup
Clone the repository using the following command:

```bash
git clone https://github.com/medunigraz/PIE-Model-Experiments.git
cd PIE-Model-Experiments
```

If available, place the unzipped `setups` and `calibration` directories from [Zenodo](https://doi.org/10.5281/zenodo.17198150) into the root directory of the repository.

Note: the experimental results were computed using the PIE-Solver executable version 1.0 and requires proprietary libraries which cannot be publicly disclosed. Binary executables compiled for the major platforms (Linux, Mac OSX or Windows) can be provided upon request, or within a complete Docker image, for replicating the results of this study by contacting the corresponding authors.

## Experimental Results
Running the following command creates the `results` output directory reproduces the experimental results:

```bash
./eval.sh -np=32
```

where the option `-np=<int>` allows to specify the number of threads used for computation. Reproducing all results, except the whole-heart model reference, should require roughly 15 minutes to compute on a desktop computer. 

https://github.com/user-attachments/assets/2b095eff-7b4e-4563-9c9f-26243b9b0e34

https://github.com/user-attachments/assets/32d53dd7-4e9e-4f32-bba9-1212fbefb81b

Computing the reference data for the human whole-heart experiment requires multiple hours, even with substantial computational resources. Our results were computed on the [Austrian Scientific Cluster](https://www.vsc.ac.at/home/). 

Whole-Heart|RD|PIE
--|--|--
<video src="https://github.com/user-attachments/assets/7dd3cb17-6cd6-4294-b481-1c3a869e1a08"></video> | <video src="https://github.com/user-attachments/assets/531b9096-f8f4-4603-a5ab-9af21fb2ea4f"></video> | <video src="https://github.com/user-attachments/assets/2bf40b8d-b9c2-4ff2-b8f9-98b26a73b1ca"></video>

https://github.com/user-attachments/assets/57e08210-89df-4337-9818-85d20d2e2e83

Run the following command to clean the repository:

```bash
./clean.sh
```

## Test System Information
Testing was conducted on an Ubuntu 22.04 LTS system equipped with 32 cores of AMD Threadripper PRO 5975WX CPU and 128GB of system memory.

## References
[1] Gillette, K., Gsell, M.A.F., Prassl, A., Karabelas, E., Reiter, U., Reiter, G., Grandits, T., Stern, D., Urschler, M., Bayer, J., Augustin, C.M., Neic, A., Pock, T., Vigmond, E., Plank, G.: A Framework for the Generation of Digital Twins of Cardiac Electrophysiology from Clinical 12-leads ECGs. Med. Imag. Anal. 71, 102080 (2021)

[2] Gsell, M.A.F., Neic, A., Bishop, M.J., Gillette, K., Prassl, A.J., Augustin, C.M., Vigmond, E.J., Plank, G.: ForCEPSS—A framework for cardiac electrophysiology simulations standardization. Comput. Methods Programs Biomed. 251, 108189 (2024)

## How to Cite
When using the data provided in this repository, please cite the paper "A Physics-Informed Eikonal Model for Simulating Arrhythmias in the Human Heart in Real-Time". 

## Licence
This repository is released under the Apache Software License, Version 2.0 ("Apache 2.0").

