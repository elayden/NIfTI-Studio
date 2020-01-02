% Define NIfTI Studio toolbox base path:
nifti_studio_path = 'NIfTI-Studio'; % Insert custom path to repo here

% Full paths to grey matter (GM) and white matter (WM) images:
GM_path = fullfile(nifti_studio_path, '\data\human\mni_icbm152_gm_tal_nlin_sym_09a.nii');
WM_path = fullfile(nifti_studio_path, '\data\human\mni_icbm152_wm_tal_nlin_sym_09a.nii');

% Threshold GM and WM to remove low intensity halo which obscures gyri/sulci, etc.
GM = load_nii(GM_path); % NIfTI_tools function
GM.img(GM.img < .5) = 0;
save_nii(GM, fullfile(nifti_studio_path, '\data\human\mni_GM_thresh.nii')) % NIfTI_tools function

%% Open image with main NIfTI Studio editor
handles = nifti_studio('background', GM); %#ok % load from previously loaded image structure

% handles = nifti_studio('background', ...
%      fullfile(nifti_studio_path, '\data\human\mni_GM_thresh.nii')) % load from file path

% Draw spherical regions of interest (ROIs):
    % Select -> New Overlay...
    % Draw -> Select Draw Color -> 1
    % Draw -> Shapes -> Sphere -> 6 (mm)
    % Click on image to generate spherical ROI
    % Save drawing as overlay: Save -> Save Current Overlay -> 'spherical_rois.nii'

%% Plot 3D rendered image w/ spherical ROIs and connectivity data:
    % Note: this cnn also be accomplished within the nifti_studio GUI:
    % Display -> Orientation -> 3D Display

% Create a connectivity matrix:
connmat = zeros(6);
connmat(4, 5) = .6; % x- and y- indices of matrix denote ROI numbers
connmat(5, 6) = .4; % specify connections between ROIs 
connmat(4, 6) = -.4; % negative numbers are by default shown as blue, positive as red (other options available)

handles = nifti_studio_3D('background', GM, ...
    'ROI', fullfile(nifti_studio_path, '\data\human\spherical_rois.nii'),...
    'connmat', connmat); %#ok

% Desired output:  human_brain_3d_rois_connections.png

%% Plot mosaic of slices:
handles = nifti_studio_mosaic('background', GM, ...
    'overlay', fullfile(nifti_studio_path, '\data\human\spherical_rois.nii'),...
    'dimension',3);

% Desired output:  human_brain_mosaic_axial.png