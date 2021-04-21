------------------------------------------------------------------------------------------------------
# Network Refinement: A random walk-based framework for network denoising
------------------------------------------------------------------------------------------------------

## Overview
NR (Network Refinement) is a network denoising framework that adjusts the edge weights by applying a graph nonlinear operator based on a diffusion process on random walks. Specifically, the framework consists of two approaches named NR-F and NR-B, which improve the SNR of noisy input networks from two different perspectives: NR-F aims at enhancing signal strength, while NR-B aims at weakening noise strength. We will choose from which angle to improve the SNR of the network according to the characteristics of the network itself. We show that NR can: (1) Improve the accuracy of Community Detection; (2) Enhance the quality of Hi-C contact map; (3) Distinguish strong and weak relationships; (4) Increase the accuracy of Gene Regulatory Networks inference.

## Running NR in Matlab or R
If you want to reproduce the results of the corresponding experiment, go to the corresponding folder and run the corresponding code. For example, if you want to test that NR-F can improve the accuracy of Community Detection, follow the instructions of <strong>Application1: NR-F improves the accuracy of Community Detection</strong>.

## Simulated noisy networks
- This contains two folders, one for the simulated experiments of NR-F and one for the simulated experiments of NR-B. 
- We give the codes for data simulation, SNR calculation, AUPR and AUROC calculation for each experiment. 
- To enable you to reproduce our results exactly, we give the simulation data of our experiments in the folders. You need to unzip them before using.

## Application1: NR-F improves the accuracy of Community Detection.
1. Open the <strong>Zachary.Rproj</strong>.
2. Run the <strong>Zachary.R</strong> to reproduce our results.

## Application2: NR-F enhances the quality of H
