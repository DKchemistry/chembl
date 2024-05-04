# Strain Metric Benchmarking

Benchmarking the Utility of Strain Metrics in Structure Based Virtual Screening

## Table of Contents

- [Introduction](#introduction)
- [Methods](#methods)
  - [Torsional Database](#torsional-database)
  - [Docking Methodology](#docking-methodology)
  - [Benchmarks Used](#benchmarks-used)
- [Results](#results)
  - [Data Analysis Pipeline](#data-analysis-pipeline)
  - [Performance Metrics](#performance-metrics)
- [Further Analysis](#further-analysis)
- [Collaboration and Future Directions](#collaboration-and-future-directions)
- [References](#references)
- [Acknowledgments](#acknowledgments)

## Introduction

This benchmarking study is meant to serve as a proof of concept work for the integration of strain metrics into the [Deep Docking platform](https://github.com/jamesgleave/DD_protocol). The identification of a ligand strain method with improved performance in structure-based virtual screening metrics that is amenable to the scale of modern ultra large virtual libraries would encourage the development of a multi-task learning approach for the Deep Docking platform. 

## Methods

### Benchmarks

Two benchmarks (GPCR-Bench and LIT-PCBA) have been used as benchmarking sets for docking and strain methods. 

GPCR-Bench is a benchmarking set of 24 GPCR targets with ~200 diverse, known actives and ~10,000 DUD-E decoys. The dataset was design by scientists at Heptares and is described here: [[GPCR-Bench: A Benchmarking Set and Practitioners’ Guide for G Protein-Coupled Receptor Docking](https://pubs.acs.org/doi/10.1021/acs.jcim.5b00660)]. The original data can be found at the following GitHub repository: [GPCR-Bench](https://github.com/dahliaweiss/GPCR-Bench)

Protein grids were already prepared via Schrodinger's Protein Preparation Wizard. The ligands are protonated SMILES. Stereochemistry was enumerated using [MayaChemTools and rdkit](http://www.mayachemtools.org/docs/scripts/html/RDKitEnumerateStereoisomers.html), analagously to Deep Docking workflows. Ligands were prepared as 3D conformers using Schrodinger's LigPrep, allowing only the previously defined protomer and stereoisomer. 

LIT-PCBA is a benchmarking set of 14 proteins (with one protein available as both the `active` and `inactive` state). The amount of actives and true inactives varied per target, as decribed in their publication: [LIT-PCBA: An Unbiased Data Set for Machine Learning and Virtual Screening](https://pubs.acs.org/doi/10.1021/acs.jcim.0c00155). The original dataset can be found at their university address: [LIT-PCBA](https://drugdesign.unistra.fr/LIT-PCBA/). 

The proteins were cleaned but not prepared for Glide docking in the original dataset. An automated workflow using Schrodinger's `X-Glide` and associated utilities were used to prepare 15 docking grids, corresponding to each protein's highest resolution PDB template and its state. The resulting grids were visually assessed for validity.

The ligands were prepared for docking using the Deep Docking workflow, as no prior protonation state was assigned. 

### Torsional Ligand Strain Method

Torsional ligand strain was calculated using "torsional database" method and code described by Gu et al (Shoichet/Irwin Groups) in [Ligand Strain Energy in Large Library Docking](https://pubs.acs.org/doi/10.1021/acs.jcim.1c00368). Minor refactoring was performed to handle easier parallel script execution and file I/O. The refactored code is available upon request and will be generally made public shortly. No changes were made to the underlying strain calculation method. 

The strain was calculated from the bound conformer and only the total ligand strain was kept. Conformer's that recieved "flagged" strain values were removed from the dataset. 

### Docking Methodology

Standard Glide SP docking was performed in Schrodinger 2021-4 using the aforementioned grids and ligands. Ligands that failed to dock were removed from subsequent analysis.

## Results

### Data Analysis Pipeline

Data analysis Jupyter Notebooks for the individual protein/ligand pairs can be found for both sets in the following directories, depending on the benchmark used:

`GPCR-Bench-master/papermill/notebooks`  

`grids_lit-pcba/papermill/notebooks`

Summary statistis regarding the average performance of the metrics for both respective benchmarking dataset can be found here: 

`GPCR-Bench-master/papermill/csv/combined_data_analysis.ipynb`

`grids_lit-pcba/papermill/csv/combined_data_analysis.ipynb`

The python package [Papermill](https://papermill.readthedocs.io/en/latest/) was used to ensure reproducibility of the data analysis pipeline. 

### Performance Metrics

In addition to the LogAUC values given for a similar analysis in Gu et al, the following metrics were also calculated: 

* EF1%
* EF5%
* AUC

We also considered the use of pareto fronts as ranking logic to evaluate performance, but preliminary results were not encouraging for further analysis. 

## Further Analysis

Code review from supervisor Prof. Francesco Gentile and MSc student, Stasa Skorupan is in progress to confirm the validity of the analysis. 

Our current findings have encouraged us to pursue this method in a multi-task learning context. 

## Collaboration and Future Directions

This work was supervised by Prof. Francesco Gentile (University of Ottawa) and my PI Prof. Thomas Frimurer (University of Copenhagen). Msc student Stasa Skorupan (University of Ottawa) provided valuable feedback and code review. 

Additionally, Stasa Skorupan and I are exploring the use of OpenEye's Freeform method as another strain metric that could be used in a similar analysis - to gauge it's utility in SBVS.

## Acknowledgments

Thank you to Prof. Francesco Gentile for his supervision and guidance, partilcuarly in ideating the project direction and verifying implementation of LogAUC metrics. Thank you to Prof. Thomas Frimurer for docking pose analysis, verying enrichment metrics,and computational support. Thank you to Stasa Skorupan for her code review and feedback.