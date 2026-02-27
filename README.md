# PIE-Model-Experiments
This repository provides a guide and the necessary resources for reproducing the experimental results of the publication: [A Physics-Informed Eikonal Model for Simulating Arrhythmias in the Human Heart in Real-Time](https://www.researchsquare.com/article/rs-7741556/v1).

## Description
The `setups` directory contains the six experimental setups from the study:
* `A1_anatomical`: an anatomically induced reentry is simulated by pacing below an isthmus-shaped scar using a S1-S2 protocol.
* `A2_functional`: a functional reentry is induced by stimulating a patch of tissue that overlaps with the waveback of a preceding planar wave front.
* `A3_wholeheart`: a ventricular tachycardia is simulated in an anatomically detailed human whole heart model that incorporates a structurally accurate infarct scar.
* `B1_restitution`: the left-hand face of a strip of healthy tissue is paced with a dynamic pacing protocol to investigate APD and CV restitution properties.
* `B2_curvature`: concave and convex wave fronts are induced in healthy tissue to investigate curvature properties.
* `B3_diffusion`: an interface between healthy and borderzone tissue types is paced at the bottom face to investigate diffusion properties.

Each setup consists of mesh files and one or multiple simulation plans in `.json` format. Meshes are stored in binary CARP (`.bpts`, `.belem`, `.blon`) and Visual Toolkit (`.vtk`) formats which can be visualized in either [NumeriCor Studio](https://numericor.at/rlb/wordpress/products/#ProdStudio) or [ParaView](https://www.paraview.org/). The human whole-heart experiment was carried over from [1] and additionally contains the lead-field data required for ECG reconstruction. 

Note: The meshes for each setup, except the `A3_wholeheart`, are generated at runtime if not present. The complete dataset including the calibration data, the generated meshes, the anonymized human whole-heart model with lead-field data and rendered videos of each experiment are available via [Zenodo](https://doi.org/10.5281/zenodo.17198150) after the embargo is lifted.

The `calibration` directory contains files with captured physiological properties that were aquired through ForCEPSS [2] and are referenced by each simulation plan. These consist of ionic model state vectors, action potential (AP) shapes, action pontential duration (APD) and conduction velocity (CV) restitution curves.

The `scripts` directory contains a collection of python scripts to generate, process or visualize the experimental data:
* `meshgen.py`: performs automatic mesh generation for an experiment in case of absence.
* `visualize_calibration.py`: visualizes the calibration results.
* `generate_lf.py`: generation of lead-field data using GIZMO.
* `dat2bdat.py`: converts lead-field data from text to binary format.
* `apply_tstart_offset.py`: corrects the timing offset of LAT files caused by prepacing in CARP.
* `error_density.py`: plots LAT, CV and LRT error density maps for a given setup.
* `track_phase_singularity.py`: performs phase singularity tracking and visualization.
* `phase_singularity_error.py`: computes quantiative phase singularity error metrics.
* `frequency_maps.py`: computes and plots dominant frequency maps from nodal PSDs.
* `phasefield_correlation.py`: computes phasemaps and phase correlation coefficients over time.
* `plot_outliers.py`: visualization of LAT, CV and LRT outliers.
* `ecg_comparison.py`: visualization of computed ECGs (legacy version).
* `ecg_comparison_v2.py`: visualization of computed ECGs (current version).
* `ecg_quantitative.py`: computes quantitative ECG metrics, including PCC, DTW and R-R intervals.
* `performance_scaling.py`: conducts the performance scaling experiment using the `A2_functional` benchmark.
* `performance_eval.py`: conducts an performance evaluation of PIE simulations.

The `results` directory contains the computed experimental data, including simulation data `sim`, paraview ensight visualizations `ens`, rendered images `png` and computed ECGs `ecg`.

The `eval.sh` bash file automatically runs all operations for reproducing the experimental data.

The `clean.sh` bash file cleans generated data within this repository.

## Setup
The evaluation relies on the [openCARP](https://opencarp.org/) simulation framework which also includes [mesher](https://git.opencarp.org/openCARP/openCARP/-/tree/master/tools/mesher) and [meshtool](https://bitbucket.org/aneic/meshtool/src/master/) required for mesh generation. Two options are presented to meet the dependencies required for reproducing the experimental results:

### Docker Setup
All neccessary prequesites for reproducing the experimental data are available within the [openCARP Docker Image](https://opencarp.org/download/installation#installation-of-opencarp-). A Dockerfile is provided within this repository to build an image with all dependencies:

```bash
docker build -t pie-img .
docker run --rm -it --shm-size=512m pie-img
```

Note: shared memory was increased to run openCARP within the container.

### Manual Setup
A [python3](https://www.python.org) environment is required to reproduce the experimental data of the PIE model. For example, use [Miniconda3](https://www.anaconda.com/docs/getting-started/miniconda/install) to create a suitable environment on your machine: 

```bash
conda create --name pie-env --file requirements.txt
conda activate pie-env
```

To install the openCARP simulation framework please refer to the official [Documentation](https://opencarp.org/download/installation) as instructions are diverse between operating systems. Finally, clone the repository using the following command:

```bash
git clone https://github.com/medunigraz/PIE-Model-Experiments.git
cd PIE-Model-Experiments
```

Note: the experimental results were computed using the PIE-Solver executable version 1.0. The executable itself cannot be publicly disclosed because it relies on proprietary libraries. However, binary executables compiled for the major platforms (Linux, macOS or Windows) can be provided upon request by contacting the corresponding authors.

## Experimental Results
Running the following command creates the `results` output directory reproduces the experimental results:

```bash
./eval.sh -np=32
```

where the option `-np=<int>` allows to specify the number of threads used for computation. Reproducing all results, except the whole-heart model reference, should take roughly 15 minutes on a desktop computer. 

https://github.com/user-attachments/assets/c4f291c4-4a62-4468-a736-08f7d34a5ccb

https://github.com/user-attachments/assets/2208a2a6-dcc6-4b3e-ae64-08208a553125

Computing the reference data for the human whole-heart experiment requires multiple hours, even with substantial computational resources. Our results were computed on the [Austrian Scientific Cluster](https://www.vsc.ac.at/home/). 

Whole-Heart|RD|PIE
--|--|--
<video src="https://github.com/user-attachments/assets/e85a4e22-a305-448b-983e-8ddf0989a511"></video> | <video src="https://github.com/user-attachments/assets/e907f6d1-2e0a-44ff-a14c-c9178dbbd4f9"></video> | <video src="https://github.com/user-attachments/assets/55ace140-d441-4259-9a6f-d2a60b9a4765"></video>

https://github.com/user-attachments/assets/cc897e84-7b0b-4035-83ca-4ba1d977bad8

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
When using the data or code provided in this repository, please cite the paper [A Physics-Informed Eikonal Model for Simulating Arrhythmias in the Human Heart in Real-Time](https://www.researchsquare.com/article/rs-7741556/v1). 

## Licence
This repository is released under the Apache Software License, Version 2.0 ("Apache 2.0").

