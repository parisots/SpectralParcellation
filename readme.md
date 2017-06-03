Multi-scale Spectral Clustering Parcellation
--------------------------------------------

 This code provides a Matlab implementation of the extension of the multi-scale normalised cuts approach (Cour et al., CVPR, 2005) for group-wise connectivity-driven brain parcellation. 
 
 Please cite any of the corresponding papers if using the code: 

-Parisot, S., Arslan, S., Passerat-Palmbach, J., Wells III, W.M., Rueckert, D.: Group-wise Parcellation of the Cortex Through Multi-scale Spectral Clustering. NeuroImage (2016)

-Parisot, S., Arslan, S., Passerat-Palmbach, J., Wells III, W.M., Rueckert, D.: Tractography-Driven Groupwise Multi-scale Parcellation of the Cortex. In: Information Processing in Medical Imaging. pp. 600–612. Springer (2015)


#### INSTALLATION
  
 In order to include dependencies folders, run the following command from the main folder when starting Matlab:
```matlab
 addpath(genpath('Dependencies'))
```

 Binaries are provided, but if necessary, recompile the normalised cut mex functions (Ncut folder) by running 'CompileDir_simple.m' in the installation directory: 
```matlab
 CompileDir_simple('./Dependencies')
```

#### USAGE

The main functions for running the code are *SpectralParcellation.m* (group-wise parcellation) and *SingleParcellation.m* (single subject parcellation).
We assume that each subject has a specific identification (it could be the file name), linked to the way the data is stored. 

We provide an example (*example.m* file) using the file organisation of the [Human Connectome Project database](https://db.humanconnectome.org) and a tractography matrix obtained from [FSL's probtrackX](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDT/UserGuide#PROBTRACKX_-_probabilistic_tracking_with_crossing_fibres) (data not provided due to the large size of the matrices). <br/>
In order to run on a different dataset, the functions *ComputeCorrelationMatrix.m* and *loadSubjectData.m* have to be replaced. 

##### Data requirements: 
- A connectivity matrix per subject (tractography or fMRI times series correlations for instance)
- A triangular cortical surface mesh, with corresponding inflated sphere. (Note that the sphere is only necessary for uniform initialisation of the seeds for supervertex parcellation.)  

#### DEPENDENCIES 

This software requires several Matlab toolboxes to run. All dependencies are provided in the 'Dependencies' folder, except for the Human Connectome Project's workbench and FSL, which can be downloaded at:<br/>
http://www.humanconnectome.org/software/connectome-workbench.html <br/>
http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation <br/>
The dependencies to the workbench, FSL and cifti toolbox can be removed if no cifti files are used in the ComputeCorrelationMatrix.m and loadSubjectData.m functions. 

Provided dependencies: 

 - [Cifti toolbox](https://github.com/Washington-University/cifti-matlab) (the link points to a more recent version)
 - [Gifti toolbox](http://www.artefact.tk/software/matlab/gifti/)
 - [Multiscale normalised cuts code](http://www.timotheecour.com/software/ncut_multiscale/ncut_multiscale.html)
 - [Fast marching toolbox](http://www.mathworks.com/matlabcentral/fileexchange/6110-toolbox-fast-marching)
 - [Iso2mesh toolbox](http://iso2mesh.sourceforge.net/cgi-bin/index.cgi)

