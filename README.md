# Estimating the Rate of Exponential Growth of SARS-CoV-2 Population

## Overview

This project aims to estimate the rate of exponential growth of the European population of SARS-CoV-2 by analyzing genome sequences. The dataset includes an observed set of SARS-CoV-2 genomes and 10,000 simulated datasets generated using different growth parameters. Each genome in the datasets is represented by a series of polymorphic positions, indicating the presence (1) or absence (0) of mutations compared to the reference genome.

## Methodology

### Data Preparation: 
The observed dataset and simulated datasets are formatted into matrices where each row represents a genome and each column represents a polymorphic position.

### Statistical Analysis:
   Pairwise Differences (K Statistic): This is the average number of differences between pairs of sequences.
        
   W Statistic: Calculated as w=S/a1​, where S is the number of polymorphic positions and a1 is a summation term based on the number of sequences.
        
   Tajima's D: A measure of the difference between the mean pairwise differences (K) and the expected number of segregating sites (W).

### Normalization: 
The K, W, and D values are normalized using the mean and variance of the corresponding statistics from the simulated datasets. The same normalization is applied to the observed dataset values.

### Distance Calculation: 
Euclidean distances between the normalized observed statistics and each normalized simulated dataset are calculated to identify the closest simulated datasets to the observed data.

### Parameter Estimation: 
The growth parameters corresponding to the 500 smallest Euclidean distances are identified, and the mean and median of these parameters are calculated to estimate the growth rate of the population.

### Visualization - Result Interpretation:
The density and histogram plots provide a clear picture of the growth rate of SARS-CoV-2 in Europe. The mean growth rate is around 110, with most values clustered close to this number. However, the presence of some extreme values indicates variability, possibly due to factors like genetic differences or public health measures.

Interestingly, the growth rate doesn’t follow a typical exponential pattern. Instead, the distribution suggests that the virus's spread may be slowing, possibly due to a decreasing number of susceptible individuals as more people become infected or gain immunity. This shift from exponential to logistic growth reflects the natural progression of an epidemic as it matures, offering valuable insights into how SARS-CoV-2 is evolving in the population.

<p align="center">
  <img src="https://github.com/user-attachments/assets/3890137d-8432-44c2-872d-1b37a41a01a6" alt="Density Plot" width="400"/>
  <img src="https://github.com/user-attachments/assets/5a2434af-fb75-4e8f-8360-66fbc542f88e" alt="Histogram" width="400"/>
</p>
