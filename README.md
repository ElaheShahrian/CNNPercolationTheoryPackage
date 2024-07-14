# Creating CNN Model using by Percolation Theory
<img width="100%" loading="lazy" src="https://github.com/SamirPaulb/SamirPaulb/blob/main/assets/rainbow-superthin.webp" />

This repository contains implementation of CNN algorithms using by percolation theory for 32x32, 64x64, 128x128 grid Networks.

<p align="center">
  <img src="https://github.com/ElaheShahrian/CNNPercolationTheoryPackage/blob/master/Sample1.png" width="800">
</p>

## Instruction

1.	Install [MATLAB](https://www.mathworks.com/).
2.	Install [Eclipse](https://www.slb.com/products-and-services/delivering-digital-at-scale/software/eclipse-industry-reference-reservoir-simulator/eclipse/).
3.	Download [CNNPercolationTheoryPackage](https://github.com/ElaheShahrian/CNNPercolationTheoryPackage) .



## Usage
### First part: Generating the images of Percolation 
#### Run Train.m including:
1.	Creating Percolation maps for 7 Probabilities (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9) (P,p)
2.	Generating Percolation maps in .INC files (k, Φ)
3.	Preparing Eclipse .DATA files
4.	Coupling MATLAB & Eclipse
5.	Reading .RSM files
6.	Calculating Darcy law’s
7.	Saving .xlsx files (Darcy permeability, binary map)

### Second part: Preprocessing binary realizations
#### Run MakeData.m including:
1.	Merging all of Probabilities for each grid Networks (32x32, 64x64, 128x128)
2.	Converting .xlsx files to .mat files


### Third part: Training Deep Learning
#### Run CNN#.m including:
1.	Loading .mat files
2.	Partitioning data (Train data, Test data )
3.	Creating CNN network layers [input layer + 2 conv + relu layer + average pooling, 2 conv + relu layer, fully connected Layer]
4.	Setting options (Mini batch size, Max epochs, Validation data, …) 
5.	Training model 
6.	 Showing results (RMSE, R2)

### Fourth part: Extrapolation (Up-scaling), Interpolation (Down-scaling)
#### Run CNNUp/Down.m including:
1.	Loading .mat files
2.	Extrapolating/ Interpolating maps for each grid Networks
3.	Partitioning data (Train data, Test data )
4.	Creating CNN network layers [input layer + 2 conv + relu layer + average pooling, 2 conv + relu layer, fully connected Layer]
5.	Setting options (Mini batch size, Max epochs, Validation data, …) 
6.	Training model 
7.	 Showing results (RMSE, R2)

### Fifth part: Results
<p align="center">
  <img src="https://github.com/ElaheShahrian/CNNPercolationTheoryPackage/blob/master/Result.jpg" width="800">
</p>


### Publication