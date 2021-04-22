------------------------------------------------------------------------------------------------------
# Network Refinement: A random walk-based framework for network denoising
------------------------------------------------------------------------------------------------------

## Overview
NR (Network Refinement) is a network denoising framework that adjusts the edge weights by applying a graph nonlinear operator based on a diffusion process on random walks. Specifically, the framework consists of two approaches named NR-F and NR-B, which improve the SNR of noisy input networks from two different perspectives: NR-F aims at enhancing signal strength, while NR-B aims at weakening noise strength. We will choose from which angle to improve the SNR of the network according to the characteristics of the network itself. We show that NR can: (1) Improve the accuracy of Community Detection; (2) Enhance the quality of Hi-C contact map; (3) Distinguish strong and weak relationships; (4) Increase the accuracy of Gene Regulatory Networks inference.

## Running NR in Matlab or R
If you want to reproduce the results of the corresponding experiment, go to the corresponding folder. For example, if you want to test that NR-F can improve the accuracy of Community Detection, follow the instructions of <strong>Application1: NR-F improves the accuracy of Community Detection</strong>.

## Simulated noisy networks
- This contains two folders, one for the simulated experiments of NR-F and one for the simulated experiments of NR-B. 
- We give the codes for data simulation, SNR calculation, AUPR and AUROC calculation for each experiment. 
- To enable you to reproduce our results exactly, we give the simulation data of our experiments in the folders. You need to unzip them before using.

## Application1: NR-F improves the accuracy of Community Detection.
1. Open the <strong>Zachary.Rproj</strong>.
2. Run the <strong>Zachary.R</strong> to reproduce our results.

## Application2: NR-F enhances the quality of Hi-C contact map from the human genome
1. You can reproduce the denoising results of NR-F on the chromosome 19  by running the <strong>test.m</strong> file. 
2. The complete results are saved in the <strong>plot result in R</strong> folder and can be displayed by running <strong>plot_performance.R</strong>.
3. If you want to reproduce all the results, please go to [https://www.ncbi.nlm.nih.gov/geo/](https://www.ncbi.nlm.nih.gov/geo/) (accession number GSE63525) and download the data yourself.

## Application3: NR-B distinguishes strong collaborations of co-authorship networks
1. Open the <strong>co-author net.Rproj</strong>.
2. Run the <strong>plot_performance.R</strong> to reproduce our results.

## Application4: NR-B improves the accuracy of Gene Regulatory Networks inference 
1. Run the <strong>test.m</strong> to reproduce our results.
