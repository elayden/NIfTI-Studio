%% Sample data for this example:

% A zebra finch brain template is available here:
% https://gin.g-node.org/elayden/Zebra_Finch_Functional_Homotopy/src/master

% This is from the following publication:
% Layden, E. A., Schertz, K. E., London, S. E., & Berman, M. G. (2019). 
% Interhemispheric functional connectivity in the zebra finch brain, 
% absent the corpus callosum in normal ontogeny. NeuroImage, 195, 113-127.

% Unzip the downloaded files in a subfolder of the main NIfTI-Studio path:
% \NIfTI-Studio\data\zebra_finch\*

%% Define NIfTI Studio toolbox paths:
nifti_studio_path = 'NIfTI-Studio'; % Insert custom path to repo here

% Full paths to grey matter (GM) and white matter (WM) images:
template = fullfile(nifti_studio_path, '\data\zebra_finch\GroupwiseTemplate.nii');
rois = fullfile(nifti_studio_path, '\data\zebra_finch\SongSystemROIs.nii');
load(fullfile(nifti_studio_path, '\data\zebra_finch\male-01_day-90_corrs.mat')) % variable 'corrmat'

% Open image with main NIfTI Studio editor
nifti_studio('background', template, 'overlay', rois)

%% Plot 3D rendered image w/ ROIs and functional connectivity data:
handles = nifti_studio_3D('background', template, ...
    'background_smoothness',2, 'background_alpha', .3, ...                      % background / template
    'ROI', rois, 'roi_alpha',1,'roi_smoothness',2, ...                          % ROIs
    'connmat', corrmat, 'connmat_thresh', [-.05, .1], 'edge_thickness',20, ...  % connections / edges
    'show_axes',0, 'effect','emphasis', 'brightness',.95, 'light_axis','y');    % general appearance

% Tip:  use up/down & left/right arrow keys to rotate brain

% See output:  zebra_finch_brain_3d_rois_connections.png

%% Plot mosaic of slices:
handles = nifti_studio_mosaic('background', template, ...
    'overlay', rois,...
    'dimension', 3);

% See output:  zebra_finch_brain_mosaic_coronal.png
