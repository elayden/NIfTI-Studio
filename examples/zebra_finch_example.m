% Define NIfTI Studio toolbox base path:
nifti_studio_path = 'C:\Users\PC-3PO\Desktop\NIfTI-Studio';

% Full paths to grey matter (GM) and white matter (WM) images:
template = fullfile(nifti_studio_path, '\data\zebra_finch\zebra_finch_template.nii');
rois = fullfile(nifti_studio_path, '\data\zebra_finch\song_system_rois.nii');

% Open image with main NIfTI Studio editor
nifti_studio('background', template, 'overlay', rois)


% Plot 3D rendered image w/ ROIs and connectivity data:



% Plot mosaic of slices:

