# NIfTI-Studio
A Matlab toolbox for NIfTI and Analyze (img/hdr) image visualization, editing, and 3D rendering

[![View NIfTI Studio on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/74729-nifti-studio)

<b>Cite:</b>

Layden, E. A. (2020). NIfTI-Studio: a Matlab toolbox for NIfTI and Analyze image visualization, editing, and 3D rendering. Zenodo. http://doi.org/10.5281/zenodo.3725006

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3725006.svg)](https://doi.org/10.5281/zenodo.3725006)

##
If you find NIfTI-Studio useful and would like to support its continued development, feel free to send a cup of coffee! :) <br><br>
[![paypal](https://www.paypalobjects.com/en_US/i/btn/btn_donateCC_LG.gif)](https://paypal.me/ElliotLayden?locale.x=en_US)

## Installation / Startup
1. Clone or download the repository
2. Open Matlab and add the NIfTI-Studio directory to Matlab's path
3. Enter nifti_studio in the command-line and press enter, or right-click nifti_studio.m -> Run
4. Select a NIfTI (.nii, .nii.gz) or Analyze (.img) image to load using the NIfTI-Studio file browser

## Usage
Github Wiki user manual coming soon!

## Enhancements and bug reports
Please create an issue and mark with the appropriate tags. If an issue concerns the main toolbox (nifti_studio.m), the "Main" tag should be used. If an issue concerns the 3D module (nifti_studio_3D.m), the "3D" tag should be used. If the issue concerns the mosaic module (nifti_studio_mosaic.m, currently in Beta), the "Mosaic" tag should be used. Additionally, labels "bug" or "enhancement" should be added. 

## Examples
### 3D rendering of human MNI brain template with regions of interest (ROIs) and connections
Code:  [human_example.m](https://github.com/elayden/NIfTI-Studio/blob/master/examples/human_example.m)
<p align="middle">
  <img align="middle" src="https://github.com/elayden/NIfTI-Studio/blob/master/examples/human_brain_3d_rois_connections.png" width="750 hspace="30" /> 
</p>    
                                                                                                                                                   
### Overlay a statistical map as colors on ROIs
Code:  [human_example.m](https://github.com/elayden/NIfTI-Studio/blob/master/examples/human_example.m)
<p align="middle">
  <img align="middle" src="https://github.com/elayden/NIfTI-Studio/blob/master/examples/human_brain_3d_rois_connections_colors.png" width="750 hspace="30" /> 
</p>                                                                                                                                            

### Slice mosaic of human MNI brain template segmented grey matter (GM) with ROIs
Code:  [human_example.m](https://github.com/elayden/NIfTI-Studio/blob/master/examples/human_example.m)
![Human Example - Mosaic](https://github.com/elayden/NIfTI-Studio/blob/master/examples/human_brain_mosaic_axial.png)

### 3D rendering of zebra finch brain template with ROIs and functional connectivity
Code:  [zebra_finch_example.m](https://github.com/elayden/NIfTI-Studio/blob/master/examples/zebra_finch_example.m)

<p float="middle">
  <img src="https://github.com/elayden/NIfTI-Studio/blob/master/examples/zebra_finch_brain_3d_rois_connections.png" width="400" />
  <img src="https://github.com/elayden/NIfTI-Studio/blob/master/examples/zebra_finch_brain_3d_rois_connections_2.png" width="400" /> 
</p>   

### Slice mosaic of zebra finch brain template with ROIs
Code:  [zebra_finch_example.m](https://github.com/elayden/NIfTI-Studio/blob/master/examples/zebra_finch_example.m)
![Zebra Finch Example - Mosaic](https://github.com/elayden/NIfTI-Studio/blob/master/examples/zebra_finch_brain_mosaic_coronal.png)
