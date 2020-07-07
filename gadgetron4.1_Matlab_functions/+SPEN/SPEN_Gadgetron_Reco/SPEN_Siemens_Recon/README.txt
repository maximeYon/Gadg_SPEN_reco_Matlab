
SPEN MRI RECONSTRUCTION PIPELINE FOR SIEMENS SCANNER (PRISMA) DATA.
 

----------------------------------------------------------------------------------------------------------

THIS SCRIPT SHOULD BE USED FOR NON-PROFIT PURPOSES ONLY. USE AT YOUR OWN RISK. CONTACT Mārtiņš Otikovs <martins.otikovs_at_weizmann.ac.il> OR Samuel Cousin <samuelfcousin_at_gmail.com> FOR ADDITIONAL INFORMATION

----------------------------------------------------------------------------------------------------------

INTRODUCTION

Ultrafast MRI is a key component in a wide range of imaging applications including functional imaging, Diffusion Weighted and Diffusion Tensor Imaging (DWI & DTI). Most of these application are based on Echo Planar Imaging (EPI); however, EPI faces known challenges when targeting heterogenous tissues, when operating at high magnetic fields, or in other instances where field inhomogeneities become important. Recently, another acquisition scheme based on Spatiotemporal Encoding (SPEN) has been proposed to collect these MRI data. A number of features make SPEN a robust alternative to EPI: it can be implemented in a fully T2*-refocused manner where inhomogeneities are compensated throughout the acquisition, its bandwidth along the blipped dimension –the more artifact prone in EPI– can be set at arbitrary values, zooming along all spatial domains is trivial, and —since SPEN’s low-bandwidth dimension is recorded directly in spatial space— each signal gives a low-resolution image that can be interleaved in a referenceless fashion to overcome motions or other artifacts. On the other hand SPEN’s processing is, unlike EPI’s, not based on a Fourier transform. This software package provides the recon pipeline needed to reconstruct multi-receive SPEN images, collected on a Siemens scanner.  
----------------------------------------------------------------------------------------------------------

GENERAL NOTES

This pipeline is written in Matlab2016a (running under Linux) and is not yet tested for other Matlab / OS versions. 
Pipeline is intended for immediate processing of data acquired on Siemens Prisma 3T scanner, but should be applicable to data acquired on any other scanner with appropriate modifications of the pipeline's data reading and ordering blocks.

----------------------------------------------------------------------------------------------------------

INSTALLATION AND DEPENDENCIES

1. This code comes together with mapVBVD for reading Siemens data and SPIRiT V0.3 that includes necessary ESPIRiT tools (http://people.eecs.berkeley.edu/~mlustig/Software.html).

2. To run the pipeline, open Main_file (located in the same folder as this README file), edit system-related paths, set location of your data and run the script.

----------------------------------------------------------------------------------------------------------

DOCUMENTATION

Hopefully you can understand what's going on by examining corresponding functions and comments accompanying them. The purpose of each processing step is described along the main function (SPEN_Reconstruction) and it's subfunctions. See Cousin et al, ISMRM 2018, #385, and Cousin et al, A regularized reconstruction pipeline for high‐definition diffusion MRI in challenging regions incorporating a per‐shot image correction, MRM, 2019, for further details.


