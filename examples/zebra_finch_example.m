% Define NIfTI Studio toolbox base path:
nifti_studio_path = 'NIfTI-Studio'; % Insert custom path to repo here

% Full paths to grey matter (GM) and white matter (WM) images:
template = fullfile(nifti_studio_path, '\data\zebra_finch\zebra_finch_template.nii');
rois = fullfile(nifti_studio_path, '\data\zebra_finch\song_system_rois.nii');
load(fullfile(nifti_studio_path, '\data\zebra_finch\male-01_day-90_corrs.mat')) % variable 'corrmat'

% Open image with main NIfTI Studio editor
nifti_studio('background', template, 'overlay', rois)

%% Plot 3D rendered image w/ ROIs and functional connectivity data:
handles = nifti_studio_3D('background', template, ...
    'background_smoothness',2, 'background_alpha', .3, ...                      % background / template
    'ROI', rois, 'roi_alpha',1,'roi_smoothness',2, ...                          % ROIs
    'connmat', corrmat, 'connmat_thresh', [-.05, .1], 'edge_thickness',20, ...  % connections / edges
    'show_axes',0, 'effect','emphasis', 'brightness',.95, 'light_axis','y');  % general appearance

% Adjust camera positioning (this was configured via manual rotation and using the Camera Toolbar under the Tools menu)
set(handles.axes,...
    'CameraTarget', [9.4117   10.4964    5.2040],...
    'CameraPosition', [8.6014   53.6719   37.2729],...
    'CameraUpVector', [-0.0032    0.5263   -1.0204],...
    'CameraViewAngle', 9.8434)

% Desired output:  zebra_finch_brain_3d_rois_connections.png

%% Plot mosaic of slices:
handles = nifti_studio_mosaic('background', template, ...
    'overlay', rois,...
    'dimension', 3);

% Desired output:  zebra_finch_brain_mosaic_coronal.png
