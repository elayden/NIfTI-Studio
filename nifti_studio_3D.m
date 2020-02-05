function [handles] = nifti_studio_3D(varargin) 
% NIfTI Studio 3D  
%   A GUI for 3D rendering NIfTI images (file types: .nii, .nii.gz, .img/.hdr)
%   
% Author:
%   Elliot Layden, The University of Chicago, 2016-19
% 
% Cite: 
%   Layden, E. A. (2019). NIfTI Studio: A Matlab Toolbox for MRI Image 
%   Editing and Visualization.
% 
% All inputs are in the form of name-value pair arguments. They may go in
% any order, as long as each input is preceded by the correct input name in
% single quotes and a comma (e.g., 'background', image_1, 'ROI', image_2, etc...).
% 
% Main Inputs:
% 'background',     Multiple input types allowed:
%                   1. char/string full-path-to or filename of 
%                   a 3D NIfTI (.nii,.nii.gz) or .img/.hdr pair (specified 
%                   referencing only the .img file)
%                        -note: It is best practice to specify the full 
%                       path to the file, using, 
%                       e.g. "background = fullfile(<path>,<filename>);".
%                   3. timeSeries can be specified as a column vector or as 
%                   a 2D matrix, wherein each column of the matrix 
%                   represents a distinct time series (useful for either 
%                   the case of multiple ROIs or multiple scans/subjects). 
%                   4. timeSeries can be specified as an already 
%                   loaded image structure as returned by load_nii or 
%                   load_untouch_nii, or as a 4D time series image matrix.
%                       -note: this option is incompatible with multiple
%                       time series
%                   5. if timeSeries is not specified or is specified as 
%                   empty [], the program will request the user to first
%                   specify the # of scans/subjects, then the user will be
%                   prompted to select the appropriate files (either as
%                   multiple 3D images per scan/subject, or as a single 4D
%                   image per scan/subject).
% 'background_colors',  same input types as background, denoting a 3D image
%                       of size equal to 'background'; values of voxels
%                       denote voxelwise colors to be applied to brain
%                       surface; color scheme can be later changed by
%                       changing the colormap or color axis in GUI
% 
% 'ROI',            Character/String: specified as a filename for a 3D 
%                   image of ROIs, denoted by sequential integer voxel 
%                   values. -Cell Array: wherein each cell contains the
%                   subscript indices of an ROI (row: voxel #; column:
%                   [x,y,z] subscripts).
% 
% 'ROI_colors',     same input types as 'ROI', denoting a 3D image
%                   of size equal to 'ROI'; values of voxels
%                   denote voxelwise colors to be applied to ROI
%                   surfaces; color scheme can be later changed by
%                   changing the colormap or color axis in GUI;
%                   'ROI_colors' should denote color values for all ROIs;
%                   input supercedes any colors denoted by
%                   'background_colors' for regions of ROIs
% 
% 'connmat',        a cell array (one cell per scan) containing a 2D
%                   matrix (rows = time points; cols = nuisance 
%                   variables); this option can be useful, e.g., if one
%                   needs to regress out subject-specific motion
%                   parameters from the data, or if one wants to
%                   investigate the effect made on correlation matrices
%                   by different preprocessing options
% 
% 'connmat_thresh', a 1x2 vector specifying (1) a threshold for negative
%                   edges (edges must be lower than this value to be
%                   displayed), and (2) a threshold for positive edges
%                   (edges must be higher than this value to be displayed)
% 
% 'titles',     
% 
% TOOLS:
% 
% Rotation: 
%   Pitch (uparrow,downarrow) 
%   Yaw (leftarrow,rightarrow)
%   Roll ('e','r')
% Zoom: 
%   -smooth zoom (click + drag left (out) or right (in)), 
%   -click zoom (left click (zoom in), shift + left click (zoom out), right
%       click (additional options)
% Pan: click and drag, or use up/down/right/left arrow keys
% 
% ADDITIONAL OPTIONS:
% 
% File Options:
% 'insert_axes' % existing axes handle(s) (must be a vector of equal length to size(connmat,3)); will plot inside these axes instead of creating new
% 'print' % char: filepath with extension to automatically print to (allowable extensions are .png, .bmp, .tiff)
% 'print_res' % integer specifying dpi (default: 300 dpi)
% 
% % View Options (views start from view([180,0]), Anterior)
% 'pitch'   % specified in degrees
% 'yaw'     % specified in degrees
% 'roll'    % specified in degrees (settings below supercede pitch, yaw, roll)
% 'view_spec' % [n_axes x 2 matrix] specifying [azimuth, elevation] (default: view(3), Matlab's default 3D view)
% 'view_angle' % [n_axes x 1 vector] see camva.m for details (can be used to zoom)
% 'cam_pos' % [n_axes x 3 (x,y,z) matrix] see campos.m for details (can be used to pan)
% 'cam_roll' % [n_axes x 1 vector] see camroll.m (can be used for rolling image, as in pitch-yaw-roll) 
% 
% % Display Options:
% 'roi_labels' % [numrois x 1 cell array] specifying roi labels (not yet enabled)
% 'roi_colors' % [numrois x 3 matrix] specifying RGB color values for each ROI
% 'roi_cmap'   % string denoting colormap over which to distribute ROI 
%               colors (if used, 'roi_colors' input is ignored; default: 
%               'jet')
% 'edge_cmap'  % string denoting colormap over which to distribute ROI 
%               colors (if used, 'edge_color' input is ignored; default: 
%               'Blue-White-Red')
% 'edge_clim',  % [1 x 2] vector denoting color limits for edges
%               default: [min, max]
% 'edge_thick_lim', % [1 x 2] vector denoting the min & max absolute values
%               for edge thickness mapping; both values must be >= 0, as
%               these values are used for mapping the thickness of both
%               positive & negative connections
% 'edge_color' % [1 x 3 rgb color vector]; sets all edge colors to single 
%               color
% 'edge_colors' % [size(connmat,1) x size(connmat,2) x 3 (rbg values)
%               matrix] of colors for edges; if used, connmat should be >0 
%               for any edges present and <=0 for any edges to be ignored
% 'vox_thresh' % numerical value indicating minimum voxel intensity for 
%               background image (default: 0)
% 'effect' % string specifying effect name (options: 'Normal','Emphasis',
%               'Sketch','Shiny','Metal','Flat')
% 'background_color' % [1x3 vector] RBG color vector
% 'lock_axes' % logical 1/0 (if a multi-axes plot, should rotation, pan, 
%               zoom, etc. be locked between axes?)
% 'show_axes' % logical 1/0 (display axes, including axes ticks, axes 
%               labels, etc.)
% 'show_grid' % logical 1/0 (display gridlines?)
% 'xlim',     % specify [vector: lower,upper] xlimits of plots in voxel units
% 'ylim',     % specify [vector: lower,upper] ylimits of plots in voxel units
% 'zlim',     % specify [vector: lower,upper] zlimits of plots in voxel units
% 'axes_color' % [1x3 vector] RBG color vector (default: [.2,.2,.2])
% 'grid_color' % [1x3 vector] RBG color vector (default: [.15,.15,.15])
% 'measure' % string: 'physical', 'voxel' (unit of measure for axes ticks)
% 'physical_units' % string denoting physical units (e.g., mm); if left 
%               blank, will attempt to determine based on NIfTI / .hdr input
% 'isovalue', % for forming isosurface, see isosurface.m function (default: .5)
% 'patch_reduce_factor' % after forming isosurface, the faces/vertices are
%                           reduced to a 'patch_reduce_factor' fraction of 
%                           the original number (i.e., .2 = 20% of original
%                           faces/vertices
% 
% % More Display Options:
% 'brightness' % numeric ranging from 0 to 1 (controls color of background lighting) (default: 1)
% 'light_axis' % string: 'x','y', or 'z' (controls which axis light originates from) (default: x)
% 'background_alpha' % numeric ranging from 0 (transparent) to 1 (opaque) (controls transparency of main image/background; default: .5)
% 'roi_alpha' % numeric ranging from 0 (transparent) to 1 (opaque) (controls transparency of ROIs; default: .6)
% 'edge_alpha' % numeric ranging from 0 (transparent) to 1 (opaque) (controls transparency of edges/connections; default: 1)
% 'background_smoothness' %  numeric ranging from 1 to 8 (controls extent to which image surface is smoothed; default: 5)
% 'roi_smoothness' % numeric ranging from 1 to 8 (default: 3)
% 'edge_thickness' % numeric ranging from 1 to 20 (default: 7)
% 'sphere_size' % numeric ranging from 0 to 15; (i.e., radius in voxels; default: 4)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_ticks = 10; % number of tick labels

% Waitbar:
h_wait = waitbar(0, 'Please wait...', 'Name', 'NIfTI Studio 3D');
try
    wbch = allchild(h_wait);
    jp = wbch(1).JavaPeer;
    jp.setIndeterminate(1);
catch
end

% Identify Function Path and Add Helper Scripts:
script_fullpath = mfilename('fullpath');
[script_path,~,~] = fileparts(script_fullpath);
addpath(genpath(script_path))

% Retrieve name-value pair inputs:
inputs = varargin;
parsed_inputs = struct('background',[],'background_colors',[],'ROI',[],...
    'ROI_colors',[],'connmat',[],'connmat_thresh',[0,0],'titles',[],...
    'insert_axes',[],'print',[],'print_res',300,'view_spec',[],...
    'view_angle',[],'cam_pos',[],'cam_roll',[],'roi_labels',[],...
    'roi_colors',[],'roi_cmap',[],'vox_thresh',0,'effect',[],...
    'background_color',[1,1,1],'lock_axes',1,'show_axes',1,'show_grid',0,...
    'axes_color',.2*ones(1,3),'grid_color',.15*ones(1,3),'unit_measure',[],...
    'physical_units','','brightness',1,'light_axis','x','background_alpha',...
    .5,'roi_alpha',.8,'edge_alpha',1,'background_smoothness',5,...
    'roi_smoothness',3,'edge_thickness',7,'edge_cmap',[],'edge_color',[],...
    'sphere_size',[],'roi_type',[],'menu_on',1,'pitch',0,'yaw',0,'roll',0,...
    'xlim',[],'ylim',[],'zlim',[],'edge_colors',[],'isovalue',.5,...
    'patch_reduce_factor', .2, 'edge_clim',[],'edge_thick_lim',[]); % defaults
poss_input = {'background','background_colors','ROI','ROI_colors',...
    'connmat','connmat_thresh','titles','insert_axes','print','print_res',...
    'view_spec','view_angle','cam_pos','cam_roll','roi_labels',...
    'roi_colors','roi_cmap','vox_thresh','effect','background_color',...
    'lock_axes','show_axes','show_grid','axis_color','grid_color',...
    'unit_measure','physical_units','brightness','light_axis',...
    'background_alpha','roi_alpha','edge_alpha','background_smoothness',...
    'roi_smoothness','edge_thickness','edge_cmap','edge_color',...
    'sphere_size','roi_type','menu_on','pitch','yaw','roll','xlim','ylim',...
    'zlim','edge_colors','isovalue','patch_reduce_factor','edge_clim',...
    'edge_thick_lim'}; % 47 possible inputs
input_ind = zeros(1,length(poss_input));
for i = 1:length(poss_input)
    j = find(strcmp(poss_input{i},inputs));
    if ~isempty(j)
        input_ind(i) = j;
        input1 = inputs{input_ind(i)+1};
        parsed_inputs.(poss_input{i}) = input1;
    end
end
roi_patch = [];

% Declare Main Variables:
background = parsed_inputs.background;
ROI = parsed_inputs.ROI;
roi_type = 1;
connmat = parsed_inputs.connmat;
connmat_thresh = parsed_inputs.connmat_thresh;
titles = parsed_inputs.titles;
vox_thresh = parsed_inputs.vox_thresh;
isovalue = parsed_inputs.isovalue;
patch_reduce_factor = parsed_inputs.patch_reduce_factor;
rotate_on = true;
zoom_on = false;
pan_on = false;

% Load Background Image
if ~isempty(background)
    switch class(background)
        case 'char'
            try 
                back_img = load_nii(background); 
            catch
                try
                    back_img = load_untouch_nii(background);
                    warning('Non-orthogonal shearing detected in affine matrix of background image. Loaded successfully without applying affine.')
                catch
                    error('Error:  failed to load background image')
                end
            end
            imdat = back_img.img; 
        case 'struct'
            try
                imdat = background.img;
                back_img = background;
            catch
                error('If input ''background'' is structure, should contain field ''.img''.')
            end
        otherwise
            if isnumeric(background)
                imdat = background;
                back_img = [];
            else
                error('Input ''background'' not recognized.')
            end 
            physical_units = 'unknown';
            pixdim = ones(1,3);
    end
else
    [background, background_path] = uigetfile({'*.nii';'*.nii.gz';'*.img'},...
        'Select a 3D image:','MultiSelect','off');
    if background_path==0 % Cancel
        disp('User cancelled action.'); return; 
    end
    if ischar(background)
        try 
            back_img = load_nii(fullfile(background_path,background)); 
        catch
            try
                back_img = (fullfile(background_path,background));
                warning('Non-orthogonal shearing detected in affine matrix of background image. Loaded successfully without applying affine.')
            catch
                error('Error:  failed to load background image')
            end
        end
        imdat = back_img.img; 
    end
end
dim = size(imdat);

% Determine Units:
if ~isempty(back_img)
    pixdim = back_img.hdr.dime.pixdim(2:4);
    try
        origin = back_img.hdr.hist.originator;
        switch bitand(back_img.hdr.dime.xyzt_units, 7) % see xform_nii.m, extra_nii_hdr.m
            case 1, physical_units = '(m)';
            case 2, physical_units = '(mm)';
            case 3, physical_units = '(microns)';
            otherwise, physical_units = '';
        end
    catch
        physical_units = '';
        origin = back_img.hdr.dime.dim(2:4)/2;
        warning('Failed to retrieve voxel units.')
    end
    if any(origin == 0)
       origin = back_img.hdr.dime.dim(2:4)/2;
    end
end

% Load Background Voxel-wise Color Data if Available
background_voxel_colors=[];
if ~isempty(parsed_inputs.background_colors)
    [background_voxel_colors,~,~] = load_ROI(parsed_inputs.background_colors,dim,'background_colors');
    background_voxel_colors(imdat<vox_thresh) = nan;
    background_cmin = min(background_voxel_colors(:));
    background_cmax = max(background_voxel_colors(:));
end

% Generate Background Patch Data:
imdat(imdat==0) = nan; phys_dim = pixdim.*dim;
% [x,y,z] = ndgrid(linspace(0,phys_dim(1),dim(1)),...
%     linspace(0,phys_dim(2),dim(2)),...
%     linspace(0,phys_dim(3),dim(3)));
[x,y,z] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
XYZ = reshape([x(:),y(:),z(:),ones(numel(x),1)],...
    size(x,1),size(x,2),size(x,3),[]);
if ~isempty(parsed_inputs.background_colors)
    main_patch = isosurface(XYZ(:,:,:,1),XYZ(:,:,:,2),XYZ(:,:,:,3), ...
        imdat>=vox_thresh, isovalue, background_voxel_colors);
else
    main_patch = isosurface(XYZ(:,:,:,1),XYZ(:,:,:,2),XYZ(:,:,:,3), ...
        imdat>=vox_thresh, isovalue);
end
main_patch = reducepatch(main_patch, patch_reduce_factor);
A = double(sparse(repmat(main_patch.faces,3,1),repmat(main_patch.faces,1,3), 1)>0);
A = sparse(1:size(A,1),1:size(A,1),1./sum(A,2))*A;
main_vertices = cell(1,8); main_vertices{1} = main_patch.vertices;
for n = 2:8
    main_vertices{n} = A*main_vertices{n-1}; 
end

% Check connmat input and generate axes positions:
if ~isempty(connmat)
    n_axes = size(connmat,3);
    handles.axes = zeros(1,n_axes);
    if n_axes==2
        n_rows = 1; n_cols = 2;
        row_l = 0; col_l = [0,.5]; row_s = 1; col_s = .499;
    elseif n_axes==3
        n_rows = 1; n_cols = 3;
        row_l = 0; col_l = [0,.333,.666]; row_s = 1; col_s = .332;
    elseif n_axes==4
        n_rows = 2; n_cols = 2;
        row_l = [.5,0]; col_l = [0,.5]; row_s = .499; col_s = .499;
    elseif n_axes<=6
        n_rows = 2; n_cols = 3;
        row_l = [.5,0]; col_l = [0,.333,.666]; row_s = .499; col_s = .332;
    elseif n_axes<=8
        n_rows = 2; n_cols = 4;
        row_l = [.5,0]; col_l = [0,.25,.5,.75]; row_s = .499; col_s = .249;
    elseif n_axes<=12
        n_rows = 3; n_cols = 4;
        row_l = [.666,.333,0]; col_l = [0,.25,.5,.75]; row_s = .332; col_s = .249;
    else
        warning('Max 12 axes (correlation matrices) are allowed at this time.')
        n_rows = 3; n_cols = 4; n_axes = 12;
        row_l = [.666,.333,0]; col_l = [0,.25,.5,.75]; row_s = .332; col_s = .249;
    end
    % Determine Positions:
    ax_pos = zeros(n_axes,4); count = 0;
    for row = 1:n_rows
        for col = 1:n_cols
            count = count + 1;
            ax_pos(count,:) = [col_l(col),row_l(row),col_s,row_s];
        end
    end     
else
    n_axes = 1; handles.axes = 12.32324;
end

% Load ROI Voxel-wise Color Data if Available
roi_voxel_colors=[];
if ~isempty(parsed_inputs.ROI_colors)
    [roi_voxel_colors,~,~] = load_ROI(parsed_inputs.ROI_colors,dim,'ROI_colors');
    roi_voxel_colors(imdat<vox_thresh) = nan;
%     roi_cmin = min(roi_voxel_colors(:));
%     roi_cmax = max(roi_voxel_colors(:));
end

% Load ROI Image & Create Patch Data:
if ~isempty(ROI)
    [roidat,numrois,roi_dim] = load_ROI(ROI,dim,'ROI');
    roidat(roidat==0) = nan; % roi_dim = size(roidat); numrois = max(roidat(:));
    % Create Patches:
    if numrois>0
        if all(roi_dim==dim)
            cmap = jet(100);
            roi_patch = struct('p',cell(1,numrois),'sphere',cell(1,numrois)); 
            roi_colors = zeros(numrois,3); roi_vertices = cell(numrois,8); 
            for i = 1:numrois
                if ~isempty(parsed_inputs.ROI_colors)
                    roi_patch(i).p = isosurface(XYZ(:,:,:,1),XYZ(:,:,:,2),XYZ(:,:,:,3), ...
                        roidat==i, isovalue, roi_voxel_colors);
                else
                    roi_patch(i).p = isosurface(XYZ(:,:,:,1),XYZ(:,:,:,2),XYZ(:,:,:,3), ...
                        roidat==i, isovalue);
                end
                A = double(sparse(repmat(roi_patch(i).p.faces,3,1),...
                    repmat(roi_patch(i).p.faces,1,3), 1)>0);
                A = sparse(1:size(A,1),1:size(A,1),1./sum(A,2))*A;
                roi_vertices{i,1} = roi_patch(i).p.vertices;
                for n = 2:8
                    roi_vertices{i,n} = A*roi_vertices{i,n-1}; 
                end
                roi_colors(i,:) = cmap(round((i/numrois)*100),:);
            end
        else
            warning('ROI dimensions do not match background dimensions.')
        end
    end
    handles.ROI_patches = zeros(n_axes,numrois); % Patch Object Handles
    n_conn = numrois*(numrois-1)/2;
    h_line = zeros(n_axes,n_conn);
else
    roi_patch = []; connmat = []; h_line = []; numrois = 0;
end

% Parse titles
if isempty(parsed_inputs.titles) 
    titles = cell(1,n_axes);
else
    if ischar(parsed_inputs.titles)
        titles = {parsed_inputs.titles};
    elseif iscell(parsed_inputs.titles)
        titles = parsed_inputs.titles;
    end
end

%% Additional Options

% Check for valid 'insert_axes':
if ~isempty(parsed_inputs.insert_axes)
    insert_axes = parsed_inputs.insert_axes;
    dont_use = 0;
    if length(insert_axes) ~= n_axes
        warning('Input ''insert_axes'' is not of equal length to size(connmat,3); generating new axes instead.')
        dont_use = 1;
    else % additional checking
        for i = 1:n_axes
            if ~isgraphics(insert_axes(i),'axes')
                dont_use = 1;
                warning('One (or more) axes supplied in input ''insert_axes'' is not a valid axes object; generating new axes instead.')
                break;
            end
        end
    end
    if dont_use
        handles.axes = repmat(91.38372, n_axes, 1);
        parsed_inputs.insert_axes = [];
    else
        handles.axes = insert_axes;
    end
end
     
background_color = parsed_inputs.background_color;
lock_axes = parsed_inputs.lock_axes;
show_axes = parsed_inputs.show_axes;
show_grid = parsed_inputs.show_grid;
axes_color = parsed_inputs.axes_color;
grid_color = parsed_inputs.grid_color;
if ~isempty(parsed_inputs.physical_units) && ischar(parsed_inputs.physical_units)
    physical_units = parsed_inputs.physical_units;
end
brightness = parsed_inputs.brightness;
background_alpha = parsed_inputs.background_alpha;
roi_alpha = parsed_inputs.roi_alpha;
edge_alpha = parsed_inputs.edge_alpha;
background_smoothness = parsed_inputs.background_smoothness;
roi_smoothness = parsed_inputs.roi_smoothness;
edge_thickness = parsed_inputs.edge_thickness;
sphere_size = parsed_inputs.sphere_size;
if isempty(sphere_size); sphere_size = 4; end

% Initialize Handles & Options:
handles = struct('figure',[],'axes',[],'titles',[],'background_patches',[],...
    'ROI_patches',[],'linkObj',[],'lighting',[]);
handles.figure = 10.28173;
handles.axes = repmat(10.28173,1,n_axes);
handles.titles = repmat(10.28173,1,n_axes);
handles.background_patches = repmat(10.28173,1,n_axes);
handles.ROI_patches = repmat(10.28173,1,n_axes);
handles.linkObj = 33.121348;
handles.lighting = repmat(10.385820,n_axes,2);

lighting_check = repmat(28.17373,1,3);
h_more_opts_fig = 82.28277;
h_colorbar_fig = 28.332737;
ROI_colorbar_fig = 37.163224;
h_edge_colorbar_fig = 28.17278;
background_colors = [1,1,1;0,0,0;.2,.6,.7;0,0,1;1,0,0;0,1,0;0,1,1;1,0,1;1,1,0];
roi_single_colors = [1,0,0;0,0,1;0,1,0;0,1,1;1,0,1;1,1,0;0,0,0];
roi_cmap_options = {'jet','hot','cool','hsv','bone','colorcube','copper',...
            'spring','summer','winter','pink','gray'};

%% Initialize Figure & Menus
if isempty(parsed_inputs.insert_axes)
    handles.figure = figure('menubar','none','color',background_color,'numbertitle',...
        'off','name','NIfTI Studio 3D','units','norm','Position',[.215,.099,.595,.802]);  % .25,.16,.51,.69
else
    handles.figure = get(parsed_inputs.insert_axes,'Parent');
    set(handles.figure,'color',background_color);
end

if parsed_inputs.menu_on
    file_menu = uimenu(handles.figure,'Label','File');
        uimenu(file_menu,'Label','Save Figure','Callback',@save_figure_callback);
        uimenu(file_menu,'Label','Print Image','Callback',{@print_callback,0});
        uimenu(file_menu,'Label','Write Background as .STL (3D Printing)','Callback',@write_stl)
    tool_menu = uimenu(handles.figure,'Label','Tools');
        h_tools(1) = uimenu(tool_menu,'Label','Rotate','Checked','on','Callback',{@change_tool,1});
        h_tools(2) = uimenu(tool_menu,'Label','Zoom','Callback',{@change_tool,2});
        h_tools(3) = uimenu(tool_menu,'Label','Pan','Callback',{@change_tool,3});
        uimenu(tool_menu,'Label','Camera Toolbar','Checked','off','Callback',@add_camera_toolbar);
    view_menu = uimenu(handles.figure,'Label','View');
        h_views(1) = uimenu(view_menu,'Label','Left (Y-Z)','Callback',{@change_view,1});
        h_views(2) = uimenu(view_menu,'Label','Right (Y-Z)','Callback',{@change_view,2});
        h_views(3) = uimenu(view_menu,'Label','Anterior (X-Z)','Callback',{@change_view,3});
        h_views(4) = uimenu(view_menu,'Label','Posterior (X-Z)','Callback',{@change_view,4});
        h_views(5) = uimenu(view_menu,'Label','Superior (X-Y)','Callback',{@change_view,5});
        h_views(6) = uimenu(view_menu,'Label','Inferior (X-Y)','Callback',{@change_view,6});
    display_menu = uimenu(handles.figure,'Label','Display');
        effects_menu = uimenu(display_menu,'Label','Effects');
            h_effects(1) = uimenu(effects_menu,'Label','Normal','Checked','on','Callback',{@change_effects,1});
            h_effects(2) = uimenu(effects_menu,'Label','Emphasis','Callback',{@change_effects,2});
            h_effects(3) = uimenu(effects_menu,'Label','Sketch','Callback',{@change_effects,3});
            h_effects(4) = uimenu(effects_menu,'Label','Shiny','Callback',{@change_effects,4});
            h_effects(5) = uimenu(effects_menu,'Label','Metal','Callback',{@change_effects,5});
            h_effects(6) = uimenu(effects_menu,'Label','Flat','Callback',{@change_effects,6});
        background_menu = uimenu(display_menu,'Label','Background Color');
            h_background(1) = uimenu(background_menu,'Label','White','Checked','on','Callback',{@change_background,1});
            h_background(2) = uimenu(background_menu,'Label','Black','Callback',{@change_background,2});
            h_background(3) = uimenu(background_menu,'Label','Aqua','Callback',{@change_background,3});
            h_background(4) = uimenu(background_menu,'Label','Blue','Callback',{@change_background,4});
            h_background(5) = uimenu(background_menu,'Label','Red','Callback',{@change_background,5});
            h_background(6) = uimenu(background_menu,'Label','Green','Callback',{@change_background,6});
            h_background(7) = uimenu(background_menu,'Label','Cyan','Callback',{@change_background,7});
            h_background(8) = uimenu(background_menu,'Label','Magenta','Callback',{@change_background,8});
            h_background(9) = uimenu(background_menu,'Label','Yellow','Callback',{@change_background,9});
            h_background(10) = uimenu(background_menu,'Label','More...','Callback',{@change_background,10});
        axes_menu = uimenu(display_menu,'Label','Axes Options');
            if n_axes>1; uimenu(axes_menu,'Label','Lock Axes','Checked','on','Callback',@change_lock_axes); end
            if parsed_inputs.show_axes
                uimenu(axes_menu,'Label','Show Axes','Checked','on','Callback',@change_show_axes);
            else
                uimenu(axes_menu,'Label','Show Axes','Checked','off','Callback',@change_show_axes);
            end
            if parsed_inputs.show_grid
                menu_grid = uimenu(axes_menu,'Label','Show Grid','Checked','on','Callback',@change_show_grid);
            else
                menu_grid = uimenu(axes_menu,'Label','Show Grid','Checked','off','Callback',@change_show_grid);
            end
            axes_limits_menu = uimenu(axes_menu,'Label','Axes Limits');
                uimenu(axes_limits_menu,'Label','X Limits','Callback',{@change_axes_limits,1});
                uimenu(axes_limits_menu,'Label','Y Limits','Callback',{@change_axes_limits,2});
                uimenu(axes_limits_menu,'Label','Z Limits','Callback',{@change_axes_limits,3});
            menu_axes_color = uimenu(axes_menu,'Label','Color');
                h_ax_color(1) = uimenu(menu_axes_color,'Label','Grey','Checked','on','Callback',{@change_axes_color,1});
                h_ax_color(2) = uimenu(menu_axes_color,'Label','Black','Checked','off','Callback',{@change_axes_color,2});
                h_ax_color(3) = uimenu(menu_axes_color,'Label','White','Checked','off','Callback',{@change_axes_color,3});
        measurements_menu = uimenu(display_menu,'Label','Measurements');
            h_units(1) = uimenu(measurements_menu,'Label','Physical Units','Checked','on','Callback',{@change_units,1});
            h_units(2) = uimenu(measurements_menu,'Label','Voxels','Checked','off','Callback',{@change_units,2});
        uimenu(display_menu,'Label','More...','Callback',@more_display_options);  
    if ~isempty(parsed_inputs.background_colors)
        background_menu = uimenu(handles.figure,'Label','Main Surface');
            uimenu(background_menu,'Label','Color Limits','Callback',@change_background_clim);
            background_color_spec_menu = uimenu(background_menu,'Label','Color Spectrum');
                background_color_spec(1) = uimenu(background_color_spec_menu,'Label','jet','Checked','on','Callback',{@background_cmap_callback,1});
                background_color_spec(2) = uimenu(background_color_spec_menu,'Label','Blue-White-Red','Callback',{@background_cmap_callback,2});
                background_color_spec(3) = uimenu(background_color_spec_menu,'Label','Red-Blue','Callback',{@background_cmap_callback,3});
                background_color_spec(4) = uimenu(background_color_spec_menu,'Label','hot','Callback',{@background_cmap_callback,4});
                background_color_spec(5) = uimenu(background_color_spec_menu,'Label','cool','Callback',{@background_cmap_callback,5});
                background_color_spec(6) = uimenu(background_color_spec_menu,'Label','hsv','Callback',{@background_cmap_callback,6});
                background_color_spec(7) = uimenu(background_color_spec_menu,'Label','bone','Callback',{@background_cmap_callback,7});
                background_color_spec(8) = uimenu(background_color_spec_menu,'Label','colorcube','Callback',{@background_cmap_callback,8});
                background_color_spec(9) = uimenu(background_color_spec_menu,'Label','copper','Callback',{@background_cmap_callback,9});
                background_color_spec(10) = uimenu(background_color_spec_menu,'Label','spring','Callback',{@background_cmap_callback,10});
                background_color_spec(11) = uimenu(background_color_spec_menu,'Label','summer','Callback',{@background_cmap_callback,11});
                background_color_spec(12) = uimenu(background_color_spec_menu,'Label','winter','Callback',{@background_cmap_callback,12});
                background_color_spec(13) = uimenu(background_color_spec_menu,'Label','pink','Callback',{@background_cmap_callback,13});
                background_color_spec(14) = uimenu(background_color_spec_menu,'Label','gray','Callback',{@background_cmap_callback,14});
    end
    if ~isempty(roi_patch)
        rois_menu = uimenu(handles.figure,'Label','ROIs');
        roi_type_menu = uimenu(rois_menu,'Label','Type');
            ROI_patches_types(1) = uimenu(roi_type_menu,'Label','Clusters','Checked','on','Callback',{@change_rois,1});
            ROI_patches_types(2) = uimenu(roi_type_menu,'Label','Spheres','Checked','off','Callback',{@change_rois,2,sphere_size});
        roi_colors_menu = uimenu(rois_menu,'Label','Color Scheme');
            roi_single_colors_menu = uimenu(roi_colors_menu,'Label','Single Color');
                ROI_patches_single_colors(1) = uimenu(roi_single_colors_menu,'Label','Red','Callback',{@change_roi_single_colors,1});
                ROI_patches_single_colors(2) = uimenu(roi_single_colors_menu,'Label','Blue','Callback',{@change_roi_single_colors,2});
                ROI_patches_single_colors(3) = uimenu(roi_single_colors_menu,'Label','Green','Callback',{@change_roi_single_colors,3});
                ROI_patches_single_colors(4) = uimenu(roi_single_colors_menu,'Label','Cyan','Callback',{@change_roi_single_colors,4});
                ROI_patches_single_colors(5) = uimenu(roi_single_colors_menu,'Label','Magenta','Callback',{@change_roi_single_colors,5});
                ROI_patches_single_colors(6) = uimenu(roi_single_colors_menu,'Label','Yellow','Callback',{@change_roi_single_colors,6});
                ROI_patches_single_colors(7) = uimenu(roi_single_colors_menu,'Label','Black','Callback',{@change_roi_single_colors,7});
                ROI_patches_single_colors(8) = uimenu(roi_single_colors_menu,'Label','More...','Callback',{@change_roi_single_colors,8});
            roi_color_spectrum_menu = uimenu(roi_colors_menu,'Label','Color Spectrum');
                ROI_patches_color_spec(1) = uimenu(roi_color_spectrum_menu,'Label','jet','Callback',{@roi_colormap_callback,1});
                ROI_patches_color_spec(2) = uimenu(roi_color_spectrum_menu,'Label','hot','Callback',{@roi_colormap_callback,2});
                ROI_patches_color_spec(3) = uimenu(roi_color_spectrum_menu,'Label','cool','Callback',{@roi_colormap_callback,3});
                ROI_patches_color_spec(4) = uimenu(roi_color_spectrum_menu,'Label','hsv','Callback',{@roi_colormap_callback,4});
                ROI_patches_color_spec(5) = uimenu(roi_color_spectrum_menu,'Label','bone','Callback',{@roi_colormap_callback,5});
                ROI_patches_color_spec(6) = uimenu(roi_color_spectrum_menu,'Label','colorcube','Callback',{@roi_colormap_callback,6});
                ROI_patches_color_spec(7) = uimenu(roi_color_spectrum_menu,'Label','copper','Callback',{@roi_colormap_callback,7});
                ROI_patches_color_spec(8) = uimenu(roi_color_spectrum_menu,'Label','spring','Callback',{@roi_colormap_callback,8});
                ROI_patches_color_spec(9) = uimenu(roi_color_spectrum_menu,'Label','summer','Callback',{@roi_colormap_callback,9});
                ROI_patches_color_spec(10) = uimenu(roi_color_spectrum_menu,'Label','winter','Callback',{@roi_colormap_callback,10});
                ROI_patches_color_spec(11) = uimenu(roi_color_spectrum_menu,'Label','pink','Callback',{@roi_colormap_callback,11});
                ROI_patches_color_spec(12) = uimenu(roi_color_spectrum_menu,'Label','gray','Callback',{@roi_colormap_callback,12});
    end
    if ~isempty(connmat)
    edges_menu = uimenu(handles.figure,'Label','Connections');
        uimenu(edges_menu,'Label','Color Limits','Callback',@change_edges_clim);
        edge_colors_menu = uimenu(edges_menu,'Label','Color Scheme');
        edges_single_colors_menu = uimenu(edge_colors_menu,'Label','Single Color');
            h_edges_single_colors(1) = uimenu(edges_single_colors_menu,'Label','Red','Callback',{@change_edges_single_colors,1});
            h_edges_single_colors(2) = uimenu(edges_single_colors_menu,'Label','Blue','Callback',{@change_edges_single_colors,2});
            h_edges_single_colors(3) = uimenu(edges_single_colors_menu,'Label','Green','Callback',{@change_edges_single_colors,3});
            h_edges_single_colors(4) = uimenu(edges_single_colors_menu,'Label','Cyan','Callback',{@change_edges_single_colors,4});
            h_edges_single_colors(5) = uimenu(edges_single_colors_menu,'Label','Magenta','Callback',{@change_edges_single_colors,5});
            h_edges_single_colors(6) = uimenu(edges_single_colors_menu,'Label','Yellow','Callback',{@change_edges_single_colors,6});
            h_edges_single_colors(7) = uimenu(edges_single_colors_menu,'Label','Black','Callback',{@change_edges_single_colors,7});
            h_edges_single_colors(8) = uimenu(edges_single_colors_menu,'Label','More...','Callback',{@change_edges_single_colors,8});
        edges_color_spec_menu = uimenu(edge_colors_menu,'Label','Color Spectrum');
            edges_color_spec(1) = uimenu(edges_color_spec_menu,'Label','Blue-White-Red','Checked','on','Callback',{@edges_colormap_callback,1});
            edges_color_spec(2) = uimenu(edges_color_spec_menu,'Label','Red-Blue','Callback',{@edges_colormap_callback,2});
            edges_color_spec(3) = uimenu(edges_color_spec_menu,'Label','jet','Callback',{@edges_colormap_callback,3});
            edges_color_spec(4) = uimenu(edges_color_spec_menu,'Label','hot','Callback',{@edges_colormap_callback,4});
            edges_color_spec(5) = uimenu(edges_color_spec_menu,'Label','cool','Callback',{@edges_colormap_callback,5});
            edges_color_spec(6) = uimenu(edges_color_spec_menu,'Label','spring','Callback',{@edges_colormap_callback,6});
            edges_color_spec(7) = uimenu(edges_color_spec_menu,'Label','summer','Callback',{@edges_colormap_callback,7});
            edges_color_spec(8) = uimenu(edges_color_spec_menu,'Label','winter','Callback',{@edges_colormap_callback,8});
            edges_color_spec(9) = uimenu(edges_color_spec_menu,'Label','bone','Callback',{@edges_colormap_callback,9});
            edges_color_spec(10) = uimenu(edges_color_spec_menu,'Label','copper','Callback',{@edges_colormap_callback,10});
            edges_color_spec(11) = uimenu(edges_color_spec_menu,'Label','pink','Callback',{@edges_colormap_callback,11});
            edges_color_spec(12) = uimenu(edges_color_spec_menu,'Label','gray','Callback',{@edges_colormap_callback,12});
        uimenu(edges_menu,'Label','Connection Thickness','Callback',@more_display_options)
    end
end

%% Plot:
handles.background_patches = repmat(38.3827,1,n_axes);
for axes_iter = 1:n_axes
    % Initialize Axes and Plot Background Patch
    if isempty(parsed_inputs.insert_axes)
        if ~isgraphics(handles.axes(axes_iter),'axes') % don't create if supplied by user
            if n_axes==1
                handles.axes(axes_iter) = axes('parent',handles.figure,'position',[0,0,1,1],'NextPlot','add');
            else
                handles.axes(axes_iter) = axes('Parent',handles.figure,'Position',ax_pos(axes_iter,:),'NextPlot','add');
                hold(handles.axes(axes_iter),'on')
            end
        else
            hold(handles.axes(axes_iter),'on')
        end
    else
        handles.axes(axes_iter) = parsed_inputs.insert_axes;
    end
    main_patch.vertices = main_vertices{background_smoothness};
    if ~isempty(parsed_inputs.background_colors)  
        handles.background_patches(axes_iter) = patch(main_patch,'parent',...
            handles.axes(axes_iter),'edgecolor','none','facevertexcdata',...
            main_patch.facevertexcdata,'facecolor','interp',...
            'AlphaDataMapping','none','FaceLighting','gouraud','FaceAlpha',...
            background_alpha);
%             isonormals(imdat,handles.background_patches(axes_iter)) 
    else
        handles.background_patches(axes_iter) = patch(main_patch,'parent',...
            handles.axes(axes_iter),'edgecolor','none','facevertexcdata',...
            repmat([1 1 1],size(main_patch.vertices,1),1),'facecolor','interp',...
            'AlphaDataMapping','none','FaceLighting','gouraud','FaceAlpha',...
            background_alpha); % 'gouraud' is for viewing curved surfaces
%             isonormals(imdat,handles.background_patches(axes_iter)) 
    end
%     axis(handles.axes(axes_iter),'equal'); % sets axes to larger than the image (whereas 'tight' fits tightly around)
    if isempty(parsed_inputs.view_spec)
        view(3); % Matlab's default 3D view angle
    end
    handles.lighting(axes_iter,1) = light('Position',[-1,0,0],'Visible','on','parent',handles.axes(axes_iter)); 
    handles.lighting(axes_iter,2) = light('Position',[1,0,0],'Visible','on','parent',handles.axes(axes_iter));
    set(handles.axes(axes_iter),'units','norm','color',background_color,'xcolor',...
        [.5,.5,.5],'ycolor',[.5,.5,.5],'zcolor',[.5,.5,.5]);
    xlabel(handles.axes(axes_iter),['X ',physical_units]); 
    ylabel(handles.axes(axes_iter),['Y ',physical_units]); 
    zlabel(handles.axes(axes_iter),['Z ',physical_units]);
    if show_grid; grid(handles.axes(axes_iter),'on'); 
    else
        grid(handles.axes(axes_iter),'off'); 
    end
    if show_axes; set(handles.axes(axes_iter),'Visible','on'); 
    else
        set(handles.axes(axes_iter),'Visible','off'); 
    end
    set(handles.figure,'WindowKeyPressFcn',@keypress_callback)
%     rotate3d(handles.axes(axes_iter),'on') % default tool

    % Plot ROI Patch:
    if ~isempty(roi_patch)
        hold on;
        for i = 1:numrois
            roi_patch(i).p.vertices = roi_vertices{i,roi_smoothness}; %#ok
            if ~isempty(parsed_inputs.ROI_colors)
                handles.ROI_patches(axes_iter,i) = patch(roi_patch(i).p,...
                    'Parent',handles.axes(axes_iter),'EdgeColor','none','FaceVertexCData',...
                    roi_patch(i).p.facevertexcdata,'FaceColor','interp',...
                    'AlphaDataMapping','none','FaceLighting','gouraud','FaceAlpha',roi_alpha);
            else
                handles.ROI_patches(axes_iter,i) = patch(roi_patch(i).p,...
                    'Parent',handles.axes(axes_iter),'EdgeColor','none','FaceVertexCData',...
                    repmat(roi_colors(i,:),size(roi_patch(i).p.vertices,1),1),'FaceColor','interp',...
                    'AlphaDataMapping','none','FaceLighting','gouraud','FaceAlpha',roi_alpha);
            end
        end
    end

    % 3D Line Plots (Edges):
    if nargin>2 && ~isempty(connmat)
        % Color Mapping:
        connmat1 = connmat(:,:,axes_iter);
        conndat = connmat1(triu(true(numrois),1));
        non_nan = ((conndat<=connmat_thresh(1))+(conndat>=connmat_thresh(2)))>0;
        conndat(~non_nan) = nan;
        m = 64; 
        % Color Mapping
        if ~isempty(parsed_inputs.edge_clim)
            cmin = parsed_inputs.edge_clim(1);
            cmax = parsed_inputs.edge_clim(2);
        else
            cmin = min(conndat); cmax = max(conndat);
        end
        conndat_norm = min(m,round((m-1)*(conndat-cmin)/(cmax-cmin))+1);
        conndat_norm(conndat_norm<=0) = 1; % assure no negative or 0 indices
        conndat_norm(~non_nan) = nan;
        rb_cmap = bluewhitered(m,cmin,cmax);
        conndat_cmap = zeros(numrois); conndat_cmap(triu(true(numrois),1)) = conndat_norm;
        conndat_cmap = conndat_cmap + conndat_cmap';
        % Thickness Mapping:
        conndat_abs = abs(conndat); 
        % Thickness Mapping:
        conndat_abs = abs(conndat); 
        if ~isempty(parsed_inputs.edge_thick_lim)
            cmin_abs = parsed_inputs.edge_thick_lim(1);
            cmax_abs = parsed_inputs.edge_thick_lim(2);
        else
            cmin_abs = min(conndat_abs); cmax_abs = max(conndat_abs);
        end
        conndat_thick_norm = min(edge_thickness,round((edge_thickness-1)*(conndat_abs-cmin_abs)/(cmax_abs-cmin_abs))+1);
        conndat_thick_norm(conndat_thick_norm<=0) = 1; % assure no negative or 0 indices
        conndat_thick_norm(~non_nan) = nan;
        conndat_tmap = zeros(numrois); conndat_tmap(triu(true(numrois),1)) = conndat_thick_norm;
        conndat_tmap = conndat_tmap + conndat_tmap';
        % Plot Lines:
        count = 0;
        for i = 1:numrois-1
            for j = i+1:numrois
                if ~isnan(conndat_cmap(i,j))
                    count = count+1;
                    coords = zeros(2,3);
                    coords(1,:) = mean(roi_patch(i).p.vertices);
                    coords(2,:) = mean(roi_patch(j).p.vertices);
%                         h_line(axes_iter,count) = plot3(ax1(axes_iter),coords(:,1),coords(:,2),coords(:,3),'Color',...
%                             rb_cmap(conndat_cmap(i,j),:),'LineWidth',conndat_tmap(i,j),'AlignVertexCenters','on');
                    h_line(axes_iter,count) = patchline(coords(:,1),coords(:,2),coords(:,3),...
                        'edgecolor',rb_cmap(conndat_cmap(i,j),:),...
                        'LineWidth',conndat_tmap(i,j),'Parent',handles.axes(axes_iter),'edgealpha',edge_alpha);
                end
            end
        end
    end
    if ~isempty(titles)
        for axes_iter3 = 1:n_axes
            handles.titles(axes_iter3) = title(titles{axes_iter3},'Parent',handles.axes(axes_iter3)); 
        end
    end
    set(handles.axes(axes_iter),'YDir','normal') 
end
% Lock Multiple Axes if Requested
if lock_axes && n_axes>1
    handles.linkObj = linkprop(handles.axes,{'CameraPosition','CameraTarget',...
        'CameraViewAngle','View','PlotBoxAspectRatio','XLim','YLim',...
        'ZLim','XTick','XTickLabel','YTick','YTickLabel','ZTick',...
        'ZTickLabel'}); 
end

%% Do additional things if requested by inputs:

% Change view options based on user input:
view([180,0])
if parsed_inputs.pitch~=0
    nifti_camorbit(handles.axes,parsed_inputs.pitch,0,'direction','x')  
end
if parsed_inputs.yaw~=0
    nifti_camorbit(handles.axes,parsed_inputs.yaw,0,'direction','z')  
end
if parsed_inputs.roll~=0
    nifti_camorbit(handles.axes,parsed_inputs.roll,0,'direction','y')  
end

% View Specification:
if ~isempty(parsed_inputs.view_spec) && all(size(parsed_inputs.view_spec)==[n_axes,2])
    for i = 1:n_axes
        view(handles.axes(i),parsed_inputs.view_spec(i,:))
    end
end
% View Angle
if ~isempty(parsed_inputs.view_angle) && (length(parsed_inputs.view_angle)==n_axes)
    for i = 1:n_axes
        camva(handles.axes(i),parsed_inputs.view_angle(i))
    end
end
% Camera Position
if ~isempty(parsed_inputs.cam_pos) && all(size(parsed_inputs.cam_pos)==[n_axes,3])
    for i = 1:n_axes
        campos(handles.axes(i),parsed_inputs.cam_pos(i,:))
    end
end
% Camera Roll
if ~isempty(parsed_inputs.cam_roll) && (length(parsed_inputs.cam_roll)==n_axes)
    for i = 1:n_axes
%         axes(handles.axes(i)) %#ok
%         camroll(parsed_inputs.cam_roll(i))
        roll_nifti(handles.axes(i),parsed_inputs.cam_roll(i))
    end
end

% Change ROI Colors:
if ~isempty(parsed_inputs.roi_colors) && all(size(parsed_inputs.roi_colors)==[numrois,3])
    roi_colors = parsed_inputs.roi_colors;
    change_rois([],[],2,sphere_size)
end

% Change ROI Color Spec:
if ~isempty(parsed_inputs.roi_cmap) && ischar(parsed_inputs.roi_cmap)
    cind = strfind(roi_cmap_options,parsed_inputs.roi_cmap);
    which_color = [];
    for jx = 1:length(cind)
        if ~isempty(cind{jx})
            which_color = jx;
            break;
        end
    end
    if ~isempty(which_color)
        hObject = struct('Label',parsed_inputs.roi_cmap);
        roi_colormap_callback(hObject,[],which_color)
    else
        warning('No available ROI colormap found for input ''roi_cmap''.')
    end
end

% Change Edge Colors (Single Color):
if ~isempty(parsed_inputs.edge_color)
    if length(parsed_inputs.edge_color)==3
        change_edges_single_colors([],[],[],parsed_inputs.edge_color)
    else
        warning('Input ''edge_color'' should be a single 1x3 RGB vector')
    end
end

% Change Edge Color Spec
if ~isempty(parsed_inputs.edge_cmap)
    if ischar(parsed_inputs.edge_cmap)
        hobject = struct('Label',parsed_inputs.edge_cmap); 
        switch parsed_inputs.edge_cmap
            case 'Blue-White-Red', edges_colormap_callback(hobject,[],1)
            case 'Red-Blue', edges_colormap_callback(hobject,[],2)
            case 'jet', edges_colormap_callback(hobject,[],3)
            case 'hot', edges_colormap_callback(hobject,[],4)
            case 'cool', edges_colormap_callback(hobject,[],5)
            case 'spring', edges_colormap_callback(hobject,[],6)
            case 'summer', edges_colormap_callback(hobject,[],7)
            case 'winter', edges_colormap_callback(hobject,[],8)
            case 'bone', edges_colormap_callback(hobject,[],9)
            case 'copper', edges_colormap_callback(hobject,[],10)
            case 'pink', edges_colormap_callback(hobject,[],11)
            case 'gray', edges_colormap_callback(hobject,[],12)
            otherwise
                warning('Input ''edge_colormap'' not found, using default.')
        end
    else
        warning('''edge_cmap'' should be a string denoting a valid Matlab colormap or ''Blue-White-Red'' or ''Red-Blue''')
    end
end

% Change Sphere Size:
if ~isempty(parsed_inputs.sphere_size)
    hObject = struct('Value',sphere_size);
    change_sphere_size(hObject)
end

% Change ROI Type:
if ~isempty(parsed_inputs.roi_type)
    switch parsed_inputs.roi_type
        case 'cluster'
            roi_type = 1;
            change_rois([],[],roi_type,[])
        case 'sphere'
            roi_type = 2;
            change_rois([],[],roi_type,sphere_size)
    end
end

% Change Effects (can be specific to each axes)
if ~isempty(parsed_inputs.effect) && ischar(parsed_inputs.effect) && ~strcmpi('Normal',parsed_inputs.effect)
    effect_options = {'Normal','Emphasis','Sketch','Shiny','Metal','Flat'};
    which_effect = [];
    for jx = 1:length(effect_options)
        TF = strcmpi(effect_options{jx},parsed_inputs.effect); % case insensitive
        if TF
            which_effect = jx;
            break;
        end
    end
    if ~isempty(which_effect)
        change_effects([],[],which_effect)
    else
        warning('No available effect found for input ''effect''.')
    end
end

% Change unit of measure if specified:
if ~isempty(parsed_inputs.unit_measure) && ischar(parsed_inputs.unit_measure)
    if strcmpi('Voxel',parsed_inputs.unit_measure)
        change_units([],[],2)
    elseif strcmpi('Physical',parsed_inputs.unit_measure)
    else
        warning('Invalid input for ''unit_measure''; disregarded')
    end
end

% Change Light Axis:
if ~isempty(parsed_inputs.light_axis) && ischar(parsed_inputs.light_axis)
    if strcmpi('x',parsed_inputs.light_axis)
    elseif strcmpi('y',parsed_inputs.light_axis)
        change_light([],[],2)
    elseif strcmpi('z',parsed_inputs.light_axis)
        change_light([],[],3)
    else
        warning('Invalid input for ''light_axis'' (should be ''x'',''y'', or ''z'').')
    end
end

% Change Brightness:
if brightness~=1
    if parsed_inputs.brightness>=0 && parsed_inputs.brightness<=1
        hobject = struct('Value',parsed_inputs.brightness);
        change_brightness(hobject)
    else
        warning('Invalid input for ''brightness'' (should be 0 =< brightness <= 1).')
    end
end

% Change Edge Colors if specified:
if ~isempty(parsed_inputs.edge_colors) && isnumeric(parsed_inputs.edge_colors)
    % upgrade this function
    
    % Delete Lines:
    for axes_iterx = 1:n_axes
        for ix1 = 1:size(h_line,2)
            if isgraphics(h_line(axes_iterx,ix1),'patch')
                delete(h_line(axes_iterx,ix1)); 
            end
        end
    end
    
    % Plot Lines:
    count = 0;
    for ix1 = 1:numrois-1
        for jx1 = ix1+1:numrois
            if connmat(ix1,jx1)>0
                count = count+1;
                coords = zeros(2,3);
                coords(1,:) = mean(roi_patch(ix1).p.vertices);
                coords(2,:) = mean(roi_patch(jx1).p.vertices);
                h_line(axes_iterx,count) = patchline(coords(:,1),coords(:,2),coords(:,3),...
                    'edgecolor',squeeze(parsed_inputs.edge_colors(ix1,jx1,:)),...
                    'LineWidth',5,...
                    'Parent',handles.axes(axes_iterx),'edgealpha',edge_alpha);
            end
        end
    end
    
%     for axes_iterx = 1:n_axes
%         count = 0;
%         for ix1 = 1:numrois-1
%             for jx1 = ix1+1:numrois
%                 if connmat(ix1,jx1)>0 && ~any(isnan(squeeze(parsed_inputs.edge_colors(ix1,jx1,:))))
%                     count = count+1;
%                     if isgraphics(h_line(axes_iterx,count),'patch')
%                         set(h_line(axes_iterx,count),'edgecolor',squeeze(parsed_inputs.edge_colors(ix1,jx1,:)))
%                     end
%                 end
%             end
%         end
%     end
end

% Change X-Lim, Y-Lim, Z-Lim if requested:
if ~isempty(parsed_inputs.xlim)
    change_axes_limits([],[],[],parsed_inputs.xlim,[],[]);
end
if ~isempty(parsed_inputs.ylim)
    change_axes_limits([],[],[],[],parsed_inputs.ylim,[]);
end
if ~isempty(parsed_inputs.zlim)
    change_axes_limits([],[],[],[],[],parsed_inputs.zlim);
end

% Adjust tick lables (physical or voxel)
setTickLabels;

% Print: this must go last (after all other operations)
if ~isempty(parsed_inputs.print)
    print_callback([],[],1)
end

% Delete waitbar
if exist('h_wait','var') && ~isempty(h_wait) && ishandle(h_wait)
    delete(h_wait)
end

%% Callbacks:
 
% Change Tool:
function change_tool(~, ~, which_tool)
    switch which_tool
        case 1 % Rotate
%             for axes_iter1 = 1:n_axes
%                 rotate3d(handles.axes(axes_iter1),'on')
%             end
            if ~rotate_on
                pan off; pan_on = false;
                zoom off; zoom_on = false;
                set(handles.figure,'WindowKeyPressFcn',@keypress_callback)
                rotate_on = true;
                h_tools(1).Checked = 'on';
                h_tools(2).Checked = 'off';
                h_tools(3).Checked = 'off';
            else
                rotate_on = false;
                h_tools(1).Checked = 'off';
                set(handles.figure,'WindowKeyPressFcn',[])
            end
                
        case 2 % Zoom
            if ~zoom_on
                zoom on;
                zoom_on = true;
                pan_on = false;
                rotate_on = false;
                h_tools(1).Checked = 'off';
                h_tools(2).Checked = 'on';
                h_tools(3).Checked = 'off';
            else
                zoom off;
                zoom_on = false;
                h_tools(2).Checked = 'off';
            end
        case 3 % Pan;
            if ~pan_on
                pan on;
                pan_on = true;
                rotate_on = false;
                zoom_on = false;
                h_tools(1).Checked = 'off';
                h_tools(2).Checked = 'off';
                h_tools(3).Checked = 'on';
            else
                pan off;
                pan_on = false;
                h_tools(3).Checked = 'off';
            end
    end
end

function add_camera_toolbar(hObject,~,~)
    if strcmp(hObject.Checked,'off')
        cameratoolbar('Show')
        hObject.Checked = 'on';
    else
        cameratoolbar('Close')
        hObject.Checked = 'off';
    end
end

% Save Figure Function:
function save_figure_callback(~,~,~)
    % Select filename & folder for output:
    [title1,path1] = uiputfile('*.fig','Specify filename:');
    if isnumeric(path1) && path1==0
        disp('User cancelled printing.'); return;
    end
    figure_title = fullfile(path1,title1);
    savefig(handles.figure,figure_title)
end

% Print Function:
function print_callback(~,~,manually)
    % Identify File Extension:
    ext_opts = {'*.png';'*.tiff';'*.bmp'};
    if ~manually
        [title1,path1] = uiputfile(ext_opts,'Specify filename:');
        figure_title = fullfile(path1,title1);
        if isnumeric(path1) && path1==0
            disp('User cancelled printing.'); return;
        end
    else
        figure_title = parsed_inputs.print;
    end
    [~,~,file_ext] = fileparts(figure_title);
    % Identify File Extension:
    exts = {'-dpng','-dtiff','-dbmp'};
    ext_ind = [];
    for ixxx = 1:3
        if ~isempty(strfind(exts{ixxx},file_ext(2:end)))
            ext_ind = ixxx; break;
        end
    end
    if isempty(ext_ind); ext_ind = 1; end
    % Print Figure:
    set(handles.figure,'InvertHardcopy','off','PaperPositionMode','auto')
    print(handles.figure,exts{ext_ind},['-r',num2str(parsed_inputs.print_res)],'-loose',figure_title)
    disp(['Printed figure: ',figure_title])
end

% Write Stereolithography (.STL) file Callback:
function write_stl(~,~,~)
    % Select filename & folder for output:
    [title1,path1] = uiputfile('*.stl','Specify filename:');
    if isnumeric(path1) && path1==0
        disp('User cancelled'); return;
    end
    figure_title = fullfile(path1,title1);
    stlwrite(figure_title, main_patch)
    % Convert to use Matlab's new default function (same name), but must
    % perform triangulation first using delaunayTriangulation.m or
    % triangulation.m
end

function change_view(~, ~, which_view)
    for ix = 1:6; h_views(ix).Checked = 'off'; end
    h_views(which_view).Checked = 'on';
    switch which_view
        case 1, view([-90,0]) % Left
        case 2, view([90,0]) % Right
        case 3, view([180,0]) % Anterior 
        case 4, view([0,0]) % Posterior
        case 5, view([0,90]) % Superior
        case 6, view([180,-90]) % Inferior
    end
end

function change_effects(~, ~, which_effect)
    for ix = 1:6, h_effects(ix).Checked = 'off'; end
    h_effects(which_effect).Checked = 'on';
    for axes_iter1 = 1:n_axes
        axes(handles.axes(axes_iter1)) %#ok
        switch which_effect
            case 1 % normal
                material default
                lighting gouraud
            case 2 % emphasis
                material([.1,.75,.5,1,.5])
                lighting gouraud
            case 3 % Sketch
                material([.1,1,1,.25,0])
                lighting gouraud
            case 4 % Shiny
                material([.3,.6,.9,20,1])
                lighting gouraud
            case 5 % Metal
    %             material([.3 .3 1 25 .5])
                material metal
                lighting gouraud
            case 6 % Flat
                lighting flat
                material dull
        end
    end
end

function change_background(~, ~, which_color)
    for ix = 1:10; h_background(ix).Checked = 'off'; end
    h_background(which_color).Checked = 'on';
    if which_color<10
        handles.figure.Color = background_colors(which_color,:);
        for axes_iter1 = 1:n_axes
            set(handles.axes(axes_iter1),'Color',background_colors(which_color,:))
        end
        if which_color==3
            change_axes_color([],[],2)
        end
    else
        if ~ishandle(h_colorbar_fig)
            h_colorbar_fig = figure('menu','none','Name','Choose Color','NumberTitle','off');
            h_colorbar_fig.Position(4) = 70; 
            ax_colorbar = gca; set(ax_colorbar,'Position',[0,0,1,1],'Visible','off','nextplot','add');
            c_map = jet(600); c_image = zeros(100,600,3);
            c_image(:,:,1) = repmat(c_map(:,1)',100,1);
            c_image(:,:,2) = repmat(c_map(:,2)',100,1);
            c_image(:,:,3) = repmat(c_map(:,3)',100,1);
            image(c_image,'Parent',ax_colorbar,'ButtonDownFcn',{@select_color,c_map});
            set(ax_colorbar,'XLim',[0,600],'YLim',[0,100])
        else
            figure(h_colorbar_fig)
        end
    end
end

function select_color(~,~,c_map)
    pt = get(gca,'CurrentPoint');
    xpos = round(pt(1,1));
    use_color = c_map(xpos,:);
    handles.figure.Color = use_color;
    for axes_iter1 = 1:n_axes
        set(handles.axes(axes_iter1),'Color',use_color)
    end
end

function change_lock_axes(hObject,~,~)
    if strcmp(hObject.Checked,'on')
        hObject.Checked = 'off';
        link_list = {'CameraPosition','CameraTarget',...
            'CameraViewAngle','View','PlotBoxAspectRatio','XLim','YLim',...
            'ZLim','XTick','XTickLabel','YTick','YTickLabel','ZTick',...
            'ZTickLabel'};
        for ix = 1:length(link_list)
            removeprop(handles.linkObj,link_list{ix}); 
        end
        for ix = 1:length(handles.axes)
            removetarget(handles.linkObj,handles.axes(ix));
        end
    elseif strcmp(hObject.Checked,'off')
        hObject.Checked = 'on';
        handles.linkObj = linkprop(handles.axes,{'CameraPosition','CameraTarget',...
            'CameraViewAngle','View','PlotBoxAspectRatio','XLim','YLim',...
            'ZLim','XTick','XTickLabel','YTick','YTickLabel','ZTick',...
            'ZTickLabel'}); 
    end
end

function change_show_axes(hObject,~,~)
    switch hObject.Checked
        case 'on'
            hObject.Checked = 'off'; menu_grid.Checked = 'off'; show_axes = 0; 
            for axes_iter1 = 1:n_axes
                set(handles.axes(axes_iter1),'Visible','off')
            end
        case 'off'
            hObject.Checked = 'on'; show_axes = 1; 
            for axes_iter1 = 1:n_axes
                set(handles.axes(axes_iter1),'Visible','on')
            end
            if show_grid; menu_grid.Checked = 'on'; end
    end
end

function change_show_grid(hObject,~,~)
    switch hObject.Checked
        case 'on'
            hObject.Checked = 'off'; show_grid = 0;
            for axes_iter1 = 1:n_axes
                grid(handles.axes(axes_iter1),'off')
            end
        case 'off'
            hObject.Checked = 'on'; show_grid = 1;
            for axes_iter1 = 1:n_axes
                grid(handles.axes(axes_iter1),'on')
            end
    end
end

function change_axes_color(~,~,which_color)
    for ix = 1:3; h_ax_color(ix).Checked = 'off'; end
    h_ax_color(which_color).Checked = 'on';
    switch which_color
        case 1 % grey
            grid_color = [.15,.15,.15];
            axes_color = [.2,.2,.2];
        case 2 % black
            grid_color = zeros(1,3);
            axes_color = zeros(1,3);
        case 3 % white
            grid_color = ones(1,3);
            axes_color = ones(1,3);
    end
    for axes_iter2 = 1:n_axes
        set(handles.axes(axes_iter2),'GridColor',grid_color,'MinorGridColor',grid_color,...
            'GridColorMode','manual','XColor',axes_color,'YColor',axes_color,'ZColor',axes_color)
    end
end

function setTickLabels
    xlimits = xlim(handles.axes(1));
    ylimits = ylim(handles.axes(1));
    zlimits = zlim(handles.axes(1));
    xticks = linspace(xlimits(1), xlimits(2), n_ticks);
    yticks = linspace(ylimits(1), ylimits(2), n_ticks);
    zticks = linspace(zlimits(1), zlimits(2), n_ticks);
    if strcmp(h_units(1).Checked,'on') % physical units (round to 3 signif digits)
        xtick_labs = round( ( xticks - origin(1) ) * pixdim(1), 3, 'significant');
        yticks_labs = round( ( yticks - origin(2) ) * pixdim(2), 3, 'significant');
        zticks_labs = round( ( zticks - origin(3) ) * pixdim(3), 3, 'significant');
        xlab = ['X ',physical_units];
        ylab = ['Y ',physical_units];
        zlab = ['Z ',physical_units];
    else % voxel units (round to whole numbers)
        xticks = round(xticks);
        yticks = round(yticks);
        zticks = round(zticks);
        xtick_labs = xticks;
        yticks_labs = yticks;
        zticks_labs = zticks;
        xlab = 'X'; ylab = 'Y'; zlab = 'Z';
    end
    for axes_iter1 = 1:n_axes
        set(handles.axes(axes_iter1), 'XTick', xticks, ...
            'XTickLabel',xtick_labs);
        set(handles.axes(axes_iter1), 'YTick', yticks, ...
            'YTickLabel',yticks_labs);
        set(handles.axes(axes_iter1), 'ZTick', zticks, ...
            'ZTickLabel',zticks_labs);
        xlabel(handles.axes(axes_iter1), xlab); 
        ylabel(handles.axes(axes_iter1), ylab); 
        zlabel(handles.axes(axes_iter1), zlab);
    end
end

function change_units(~,~,unit_type)
    if unit_type == 1 && strcmp(h_units(1).Checked,'off')
        set(h_units(1), 'Checked', 'on')
        set(h_units(2), 'Checked', 'off')
    elseif unit_type == 2 && strcmp(h_units(2).Checked,'off')
        set(h_units(1), 'Checked', 'off')
        set(h_units(2), 'Checked', 'on')
    end
%     if strcmp(h_units(unit_type).Checked,'off')
%         switch unit_type
%             case 1 % Physical
%                 h_units(1).Checked = 'on';
%                 h_units(2).Checked = 'off';
%                 for axes_iter1 = 1:n_axes
%                     set(handles.axes(axes_iter1),'XTickLabelMode','auto','XTickMode','auto',...
%                         'YTickLabelMode','auto','YTickMode','auto',...
%                         'ZTickLabelMode','auto','ZTickMode','auto');
%                     xlabel(handles.axes(axes_iter1),['X ',physical_units]); 
%                     ylabel(handles.axes(axes_iter1),['Y ',physical_units]); 
%                     zlabel(handles.axes(axes_iter1),['Z ',physical_units]);
%                 end
%             case 2 % Voxel
%                 h_units(1).Checked = 'off';
%                 h_units(2).Checked = 'on';
%                 for axes_iter1 = 1:n_axes
%                     xtick = round(get(handles.axes(axes_iter1),'XTick')'./pixdim(1));
%                     ytick = round(get(handles.axes(axes_iter1),'YTick')'./pixdim(2));
%                     ztick = round(get(handles.axes(axes_iter1),'ZTick')'./pixdim(3));
%                     set(handles.axes(axes_iter1),'XTickLabel',num2cell(xtick),'XTickMode','manual');
%                     set(handles.axes(axes_iter1),'YTickLabel',num2cell(ytick),'YTickMode','manual');
%                     set(handles.axes(axes_iter1),'ZTickLabel',num2cell(ztick),'ZTickMode','manual');
%                     xlabel(handles.axes(axes_iter1),'X'); 
%                     ylabel(handles.axes(axes_iter1),'Y'); 
%                     zlabel(handles.axes(axes_iter1),'Z');
%                 end
%         end
%     end
    
    setTickLabels;
end

% Change Limits:
function change_axes_limits(~,~,which_axis,xlim,ylim,zlim)
    if nargin<4 % Menu Callback Functionality
        switch which_axis
            case 1 % X-axis
                prompt = {sprintf('Specify X-Limits in voxels \n\n Lower:'),'Upper:'};
                dlg_title = 'X-Limits'; num_lines = [1,17;1,17]; 
                answer = inputdlg(prompt,dlg_title,num_lines);
                if isempty(answer); return; end
                xlim = [str2double(answer{1}),str2double(answer{2})];
                for ix = 1:n_axes
                    set(handles.axes(ix),'XLim',xlim);
                end
            case 2 % Y-axis
                prompt = {sprintf('Specify Y-Limits in voxels \n\n Lower:'),'Upper:'};
                dlg_title = 'Y-Limits'; num_lines = [1,17;1,17]; 
                answer = inputdlg(prompt,dlg_title,num_lines);
                if isempty(answer); return; end
                ylim = [str2double(answer{1}),str2double(answer{2})];
                for ix = 1:n_axes
                    set(handles.axes(ix),'YLim',ylim);
                end 
            case 3 % Z-axis
                prompt = {sprintf('Specify Z-Limits in voxels \n\n Lower:'),'Upper:'};
                dlg_title = 'Z-Limits'; num_lines = [1,17;1,17]; 
                answer = inputdlg(prompt,dlg_title,num_lines);
                if isempty(answer); return; end
                zlim = [str2double(answer{1}),str2double(answer{2})];
                for ix = 1:n_axes
                    set(handles.axes(ix),'ZLim',zlim);
                end
        end
    else % Command-line Name-Value Pair functionality
        if ~isempty(xlim)            
            if ~length(xlim)==2
                error('Input ''xlim'' should be a vector of length 2.')
            end
            for ix = 1:n_axes
                set(handles.axes(ix),'XLim',xlim); % .*pixdim(1)
            end
        end
        if ~isempty(ylim)
            if ~length(ylim)==2
                error('Input ''ylim'' should be a vector of length 2.')
            end
            for ix = 1:n_axes
                set(handles.axes(ix),'YLim',ylim); % .*pixdim(2)
            end 
        end 
        if ~isempty(zlim)
            if ~length(zlim)==2
                error('Input ''zlim'' should be a vector of length 2.')
            end
            for ix = 1:n_axes
                set(handles.axes(ix),'ZLim',zlim); % .*pixdim(3)
            end
        end 
    end
    
    setTickLabels;
end

function more_display_options(~,~,~)
    if ~isgraphics(h_more_opts_fig, 'Figure')
        h_more_opts_fig = figure('menu','none','Name','Display Options','NumberTitle',...
            'off','Color',[1,1,1],'units','norm','Position',[.03,.16,.21,.725]);
        set(gca,'Visible','off')
        uicontrol(h_more_opts_fig,'Style','text','String',...
            'Brightness','Units','normalized','BackgroundColor',[1,1,1],...
            'FontSize',10,'FontWeight','bold','Position',[.05,.94,.9,.05]);
        uicontrol(h_more_opts_fig,'Style','slider','Units',...
            'normalized','Position',[.05,.9,.9,.05],'Min',0,'Max',1,'Value',brightness,...
            'BackgroundColor',[1,1,1],'Callback',@change_brightness);
        uicontrol(h_more_opts_fig,'Style','text','String',...
            'Light Axis','Units','normalized','BackgroundColor',[1,1,1],...
            'FontSize',10,'FontWeight','bold','Position',[.05,.84,.9,.05]);
        xpos = linspace(.1,.85,3);
        light_bg = uibuttongroup('Visible','on','Units','normalized','Position',...
            [.05,.76,.9,.1],'BackgroundColor',ones(1,3),'BorderType','none','SelectionChangedFcn',@change_light);
        lighting_check(1) = uicontrol(light_bg,'Style','radiobutton','Units',...
            'normalized','Position',[xpos(1),.25,.06,.6],'HorizontalAlignment','center',...
            'BackgroundColor',[1,1,1],'String','');
        lighting_check(2) = uicontrol(light_bg,'Style','radiobutton','Units',...
            'normalized','Position',[xpos(2),.25,.06,.6],'HorizontalAlignment','center',...
            'BackgroundColor',[1,1,1],'String','');
        lighting_check(3) = uicontrol(light_bg,'Style','radiobutton','Units',...
            'normalized','Position',[xpos(3),.25,.06,.6],'HorizontalAlignment','center',...
            'BackgroundColor',[1,1,1],'String','');
        switch parsed_inputs.light_axis
            case 'x'
                set(light_bg,'SelectedObject',lighting_check(1))
            case 'y'
                set(light_bg,'SelectedObject',lighting_check(2))
            case 'z'
                set(light_bg,'SelectedObject',lighting_check(3))
        end
        uicontrol(h_more_opts_fig,'Style','text','String',...
            'X','Units','normalized','BackgroundColor',[1,1,1],'HorizontalAlignment','center',...
            'FontSize',8,'FontWeight','bold','Position',[xpos(1)+.0495,.83,.03,.025]);
        uicontrol(h_more_opts_fig,'Style','text','String',...
            'Y','Units','normalized','BackgroundColor',[1,1,1],'HorizontalAlignment','center',...
            'FontSize',8,'FontWeight','bold','Position',[xpos(2)+.0112,.83,.03,.025]);
        uicontrol(h_more_opts_fig,'Style','text','String',...
            'Z','Units','normalized','BackgroundColor',[1,1,1],'HorizontalAlignment','center',...
            'FontSize',8,'FontWeight','bold','Position',[xpos(3)-.0235,.83,.03,.025]);
        uicontrol(h_more_opts_fig,'Style','text','String',...
            'Background Opacity','Units','normalized','BackgroundColor',[1,1,1],...
            'FontSize',10,'FontWeight','bold','Position',[.05,.76,.9,.03]);
        uicontrol(h_more_opts_fig,'Style','slider','Units',...
            'normalized','Position',[.05,.7,.9,.05],'Min',0,'Max',1,'Value',background_alpha,...
            'BackgroundColor',[1,1,1],'Callback',@change_background_alpha);
        uicontrol(h_more_opts_fig,'Style','text','String',...
            'ROI Opacity','Units','normalized','BackgroundColor',[1,1,1],...
            'FontSize',10,'FontWeight','bold','Position',[.05,.66,.9,.03]);
        uicontrol(h_more_opts_fig,'Style','slider','Units',...
            'normalized','Position',[.05,.6,.9,.05],'Min',0,'Max',1,'Value',roi_alpha,...
            'BackgroundColor',[1,1,1],'Callback',@change_roi_alpha);
        uicontrol(h_more_opts_fig,'Style','text','String',...
            'Connection Opacity','Units','normalized','BackgroundColor',[1,1,1],...
            'FontSize',10,'FontWeight','bold','Position',[.05,.56,.9,.03]);
        uicontrol(h_more_opts_fig,'Style','slider','Units',...
            'normalized','Position',[.05,.5,.9,.05],'Min',0,'Max',1,'Value',edge_alpha,...
            'BackgroundColor',[1,1,1],'Callback',@change_edge_alpha);
        uicontrol(h_more_opts_fig,'Style','text','String',...
            'Background Smoothness','Units','normalized','BackgroundColor',[1,1,1],...
            'FontSize',10,'FontWeight','bold','Position',[.05,.46,.9,.03]);
        uicontrol(h_more_opts_fig,'Style','slider','Units',...
            'normalized','Position',[.05,.4,.9,.05],'Min',1,'Max',8,'Value',background_smoothness,...
            'BackgroundColor',[1,1,1],'SliderStep',[1/7,1/7],'Callback',@change_back_smooth);
        uicontrol(h_more_opts_fig,'Style','text','String',...
            'ROI Smoothness','Units','normalized','BackgroundColor',[1,1,1],...
            'FontSize',10,'FontWeight','bold','Position',[.05,.36,.9,.03]);
        uicontrol(h_more_opts_fig,'Style','slider','Units',...
            'normalized','Position',[.05,.3,.9,.05],'Min',1,'Max',8,'Value',roi_smoothness,...
            'BackgroundColor',[1,1,1],'SliderStep',[1/7,1/7],'Callback',@change_roi_smooth);
        uicontrol(h_more_opts_fig,'Style','text','String',...
            'Connection Thickness','Units','normalized','BackgroundColor',[1,1,1],...
            'FontSize',10,'FontWeight','bold','Position',[.05,.26,.9,.03]);
        uicontrol(h_more_opts_fig,'Style','slider','Units',...
            'normalized','Position',[.05,.2,.9,.05],'Min',1,'Max',35,'Value',edge_thickness,...
            'BackgroundColor',[1,1,1],'SliderStep',[.05,.05],'Callback',@change_edge_thick);
        uicontrol(h_more_opts_fig,'Style','text','String',...
            'Sphere Size','Units','normalized','BackgroundColor',[1,1,1],...
            'FontSize',10,'FontWeight','bold','Position',[.05,.16,.9,.03]);
        uicontrol(h_more_opts_fig,'Style','slider','Units',...
            'normalized','Position',[.05,.1,.9,.05],'Min',0,'Max',15,'Value',sphere_size,...
            'BackgroundColor',[1,1,1],'SliderStep',[.025,.025],'Callback',@change_sphere_size);
    else
        figure(h_more_opts_fig)
    end
end

function change_brightness(hObject,~,~)
    brightness = hObject.Value;
    for ix = 1:n_axes
%         brighten(ax1(ix),hObject.Value) 
        set(handles.axes(ix),'AmbientLightColor',[brightness,brightness,brightness])
        for iy = 1:length(handles.lighting)
            if ishandle(handles.lighting(ix,iy))
                set(handles.lighting(ix,iy),'Color',[brightness,brightness,brightness]);
            end
        end
    end
end

function change_light(~,~,which_lighting)
    if nargin<3 || isempty(which_lighting)
        which_lighting = zeros(1,3);
        for ixy = 1:3; which_lighting(ixy) = get(lighting_check(ixy),'Value'); end
        which_lighting = find(which_lighting);
    end
    for axes_iter1 = 1:n_axes
        switch which_lighting
            case 1
                parsed_inputs.light_axis = 'x';
                set(handles.lighting(axes_iter1,1),'Position',[-1,0,0]); 
                set(handles.lighting(axes_iter1,2),'Position',[1,0,0]);
            case 2
                parsed_inputs.light_axis = 'y';
                set(handles.lighting(axes_iter1,1),'Position',[0,-1,0]); 
                set(handles.lighting(axes_iter1,2),'Position',[0,1,0]);
            case 3
                parsed_inputs.light_axis = 'z';
                set(handles.lighting(axes_iter1,1),'Position',[0,0,-1]); 
                set(handles.lighting(axes_iter1,2),'Position',[0,0,1]); 
        end
    end
end

function change_background_alpha(hObject,~,~)
    background_alpha = hObject.Value;
    for axes_iter1 = 1:n_axes
        set(handles.background_patches(axes_iter1),'FaceAlpha',background_alpha);
    end
end

function change_roi_alpha(hObject,~,~)
    roi_alpha = hObject.Value;
    for axes_iter1 = 1:n_axes
        for ix = 1:numrois
            set(handles.ROI_patches(axes_iter1,ix),'FaceAlpha',roi_alpha);
        end
    end
end

function change_back_smooth(hObject,~,~)
    background_smoothness = round(hObject.Value);
    for axes_iter1 = 1:n_axes
        main_patch.vertices = main_vertices{background_smoothness};
        delete(handles.background_patches(axes_iter1))
        if ~isempty(parsed_inputs.background_colors)
            handles.background_patches(axes_iter1) = patch(main_patch,...
                'parent',handles.axes(axes_iter1),'edgecolor','none','facevertexcdata',...
                main_patch.facevertexcdata,'facecolor','interp','AlphaDataMapping','none',...
                'FaceLighting','gouraud','FaceAlpha',background_alpha); % 'gouraud' is for viewing curved surfaces
%             isonormals(imdat,handles.background_patches(axes_iter)) 
        else
            handles.background_patches(axes_iter1) = patch(main_patch,...
                'parent',handles.axes(axes_iter1),'edgecolor','none','facevertexcdata',...
                repmat([1 1 1],size(main_patch.vertices,1),1),'facecolor','interp',...
                'AlphaDataMapping','none','FaceLighting','gouraud','FaceAlpha',background_alpha); % 'gouraud' is for viewing curved surfaces
%             isonormals(imdat,handles.background_patches(axes_iter)) 
        end
    end
end

function change_roi_smooth(hObject,~,~)
    roi_smoothness = round(hObject.Value);
    for axes_iter1 = 1:n_axes
        for ix = 1:numrois
        delete(handles.ROI_patches(axes_iter1,ix))
        roi_patch(ix).p.vertices = roi_vertices{ix,roi_smoothness};
        if ~isempty(parsed_inputs.ROI_colors)
            handles.ROI_patches(axes_iter1,ix) = patch(roi_patch(ix).p,...
                'Parent',handles.axes(axes_iter1),'EdgeColor','none','FaceVertexCData',...
                roi_patch(ix).p.facevertexcdata,'FaceColor','interp',...
                'AlphaDataMapping','none','FaceLighting','gouraud','FaceAlpha',roi_alpha);
        else
            handles.ROI_patches(axes_iter1,ix) = patch(roi_patch(ix).p,...
                'Parent',handles.axes(axes_iter1),'EdgeColor','none','FaceVertexCData',...
                repmat(roi_colors(ix,:),size(roi_patch(ix).p.vertices,1),1),'FaceColor','interp',...
                'AlphaDataMapping','none','FaceLighting','gouraud','FaceAlpha',roi_alpha);
        end
%         isonormals(roidat==axes_iter1,handles.ROI_patches(axes_iter1,ix)) 
        end
    end
end

function change_edge_thick(hObject,~,~)
    edge_thickness = hObject.Value;
    for axes_iter1 = 1:n_axes
        connmat1 = connmat(:,:,axes_iter1);
        conndat = connmat1(triu(true(numrois),1));
        non_nan = ((conndat<=connmat_thresh(1))+(conndat>=connmat_thresh(2)))>0;
        conndat(~non_nan) = nan;
        conndat_abs = abs(conndat); cmin_abs = min(conndat_abs); cmax_abs = max(conndat_abs);
        conndat_thick_norm = min(edge_thickness,round((edge_thickness-1)*(conndat_abs-cmin_abs)/(cmax_abs-cmin_abs))+1);
        conndat_thick_norm(conndat_thick_norm<=0) = 1; % assure no negative or 0 indices
        conndat_thick_norm(~non_nan) = nan;
        conndat_tmap = zeros(numrois); 
        conndat_tmap(triu(true(numrois),1)) = conndat_thick_norm;
        conndat_tmap = conndat_tmap + conndat_tmap';
        count = 0;
        for ix = 1:numrois-1
            for jxx = ix+1:numrois
                if ~isnan(conndat_tmap(ix,jxx))
                    count = count+1;
                    if isgraphics(h_line(axes_iter1,count))
                        set(h_line(axes_iter1,count),'LineWidth',conndat_tmap(ix,jxx))
                    end
                end
            end
        end
    end
end

function change_sphere_size(hObject,~,~)
    sphere_size = hObject.Value;
    change_rois([],[],2,sphere_size)
end

function change_rois(~,~,roi_type,s_radius)
    if ~isempty(roi_patch)
        switch roi_type
            case 1 % Clusters
                ROI_patches_types(1).Checked = 'on'; ROI_patches_types(2).Checked = 'off';
                for axes_iter1 = 1:n_axes
                    hold on;
                    for ix = 1:numrois
                        if isgraphics(handles.ROI_patches(axes_iter1,ix)); delete(handles.ROI_patches(axes_iter1,ix)); end
                        handles.ROI_patches(axes_iter1,ix) = patch(roi_patch(ix).p,'Parent',handles.axes(axes_iter1),'EdgeColor','none','FaceVertexCData',...
                            repmat(roi_colors(ix,:),size(roi_patch(ix).p.vertices,1),1),'FaceColor','interp',...
                            'AlphaDataMapping','none','FaceLighting','gouraud','FaceAlpha',roi_alpha);
%                         isonormals(roidat==axes_iter1,handles.ROI_patches(axes_iter1,ix)) 
                    end
                end
            case 2 % Spheres
                ROI_patches_types(1).Checked = 'off'; ROI_patches_types(2).Checked = 'on';
                for axes_iter1 = 1:n_axes
                        hold on;
                    for ix = 1:numrois
                        [X,Y,Z] = sphere(64);
                        s_center = mean(roi_patch(ix).p.vertices);
                        s_radius_phys = (s_radius./dim).*phys_dim;
                        s_radius_phys = mean(s_radius_phys); % average voxel dim radius
                        X = X.*s_radius_phys; Y = Y.*s_radius_phys; Z = Z.*s_radius_phys;
                        X = X + s_center(1); Y = Y + s_center(2); Z = Z + s_center(3);
                        if isgraphics(handles.ROI_patches(axes_iter1,ix)); delete(handles.ROI_patches(axes_iter1,ix)); end
                        handles.ROI_patches(axes_iter1,ix) = surf(X,Y,Z,'Parent',handles.axes(axes_iter1),...
                            'EdgeColor','none','FaceColor',roi_colors(ix,:),...
                            'AlphaDataMapping','none','FaceLighting','gouraud','FaceAlpha',roi_alpha);
%                         isonormals(roidat==ix,handles.ROI_patches(axes_iter1,ix)) 
                    end
                end
        end
    end
end

function change_background_clim(~,~,~)
    prompt = {sprintf('Specify color limits: \n\nLower Bound:'),'Upper Bound:'};
    dlg_title = 'CLim'; num_lines = [1,25;1,25]; 
    answer = inputdlg(prompt,dlg_title,num_lines);
    if isempty(answer); return; end
    background_cmin = str2double(answer{1});
    background_cmax = str2double(answer{2});
    caxis(handles.axes,[background_cmin,background_cmax])
end

function background_cmap_callback(hObject,~,which_spec)
    m = 64;
    switch hObject.Label
        case 'Blue-White-Red'
            background_cmap = bluewhitered(m,background_cmin,background_cmax);
        case 'Red-Blue'
            background_cmap = redblue(m);
        otherwise
            background_cmap = eval([hObject.Label,'(',num2str(m),')']);
    end
    for ix = 1:length(background_color_spec)
        background_color_spec(ix).Checked = 'off';
    end; background_color_spec(which_spec).Checked = 'on';
    colormap(handles.axes,background_cmap)
end

function change_roi_single_colors(~,~,which_color)
    for ix = 1:8; ROI_patches_single_colors(ix).Checked = 'off'; end
    ROI_patches_single_colors(which_color).Checked = 'on';
    if which_color < 8
        for axes_iter1 = 1:n_axes
            for ix = 1:size(handles.ROI_patches,2)
                if strcmp(ROI_patches_types(1).Checked,'on')
                    set(handles.ROI_patches(axes_iter1,ix),'FaceVertexCData',repmat(roi_single_colors(which_color,:),size(roi_patch(ix).p.vertices,1),1))
                else
                    set(handles.ROI_patches(axes_iter1,ix),'FaceColor',roi_single_colors(which_color,:))
                end
            end
        end
    else
        if ~ishandle(ROI_colorbar_fig)
            ROI_colorbar_fig = figure('menu','none','Name','Choose Color','NumberTitle','off');
            ROI_colorbar_fig.Position(4) = 70; 
            ax_colorbar = gca; set(ax_colorbar,'Position',[0,0,1,1],'Visible','off','nextplot','add');
            c_map = jet(600); c_image = zeros(100,600,3);
            c_image(:,:,1) = repmat(c_map(:,1)',100,1);
            c_image(:,:,2) = repmat(c_map(:,2)',100,1);
            c_image(:,:,3) = repmat(c_map(:,3)',100,1);
            image(c_image,'Parent',ax_colorbar,'ButtonDownFcn',{@manual_select_roi_color,c_map});
            set(ax_colorbar,'XLim',[0,600],'YLim',[0,100])
        else
            figure(ROI_colorbar_fig)
        end
    end
end

function manual_select_roi_color(~,~,c_map)
    pt = get(gca,'CurrentPoint');
    xpos = round(pt(1,1));
    use_color = c_map(xpos,:);
    if roi_type==1
        for axes_iter1 = 1:n_axes
            for ix = 1:size(handles.ROI_patches,2)
                set(handles.ROI_patches(axes_iter1,ix),'FaceVertexCData',repmat(use_color,size(roi_patch(ix).p.vertices,1),1))
                if strcmp(ROI_patches_types(1).Checked,'on')
                    set(handles.ROI_patches(axes_iter1,ix),'FaceVertexCData',repmat(use_color,size(roi_patch(ix).p.vertices,1),1))
                else
                    set(handles.ROI_patches(axes_iter1,ix),'FaceColor',use_color)
                end
            end
        end
    elseif roi_type==2
       for axes_iter1 = 1:n_axes
            for ix = 1:size(handles.ROI_patches,2)
                set(handles.ROI_patches(axes_iter1,ix),'FaceColor',use_color)
            end
        end
    end
end

function roi_colormap_callback(hObject,~,which_color)
    for ix = 1:12; ROI_patches_color_spec(ix).Checked = 'off'; end
    ROI_patches_color_spec(which_color).Checked = 'on';
%     cmap_name = get(hObject,'Label');
    cmap_name = hObject.Label;
    cmap = eval([cmap_name,'(100)']);
    for axes_iter1 = 1:n_axes
        for ix = 1:size(handles.ROI_patches,2)
            roi_colors(ix,:) = cmap(round((ix/numrois)*100),:);
            if strcmp(ROI_patches_types(1).Checked,'on')
                set(handles.ROI_patches(axes_iter1,ix),'FaceVertexCData',repmat(roi_colors(ix,:),size(roi_patch(ix).p.vertices,1),1))
            else
                set(handles.ROI_patches(axes_iter1,ix),'FaceColor',roi_colors(ix,:))
            end
        end
    end
end

function change_edges_clim(~,~,from_other)
    if nargin<3 || isempty(from_other)
        prompt = {sprintf('Specify color limits: \n\nLower Bound (Negative Threshold):'),'Upper Bound (Positive Threshold):'};
        dlg_title = 'CLim'; num_lines = [1,35;1,35]; 
        answer = inputdlg(prompt,dlg_title,num_lines);
        if isempty(answer); return; end
        connmat_thresh(1) = str2double(answer{1}); connmat_thresh(2) = str2double(answer{2});
    end
    % Delete Lines:
    for axes_iter1 = 1:n_axes
        for ixx = 1:size(h_line,2)
            if isgraphics(h_line(axes_iter1,ixx),'patch'); delete(h_line(axes_iter1,ixx)); end
        end
    end
    % Re-Plot:
    for axes_iter1 = 1:n_axes
        connmat1 = connmat(:,:,axes_iter1);
        conndat = connmat1(triu(true(numrois),1));
        non_nan = ((conndat<=connmat_thresh(1))+(conndat>=connmat_thresh(2)))>0;
        conndat(~non_nan) = nan;
        % Color Mapping
        if ~isempty(parsed_inputs.edge_clim)
            cmin = parsed_inputs.edge_clim(1);
            cmax = parsed_inputs.edge_clim(2);
        else
            cmin = min(conndat); cmax = max(conndat);
        end
        if strcmp(edges_color_spec(1).Checked,'on'); rb_cmap = bluewhitered(m,cmin,cmax); end
        conndat_norm = min(m,round((m-1)*(conndat-cmin)/(cmax-cmin))+1);
        conndat_norm(conndat_norm<=0) = 1; % assure no negative or 0 indices
        conndat_norm(~non_nan) = nan;
        conndat_cmap = zeros(numrois); conndat_cmap(triu(true(numrois),1)) = conndat_norm;
        conndat_cmap = conndat_cmap + conndat_cmap';
        % Thickness Mapping:
        conndat_abs = abs(conndat); 
        if ~isempty(parsed_inputs.edge_thick_lim)
            cmin_abs = parsed_inputs.edge_thick_lim(1);
            cmax_abs = parsed_inputs.edge_thick_lim(2);
        else
            cmin_abs = min(conndat_abs); cmax_abs = max(conndat_abs);
        end
        conndat_thick_norm = min(edge_thickness,round((edge_thickness-1)*(conndat_abs-cmin_abs)/(cmax_abs-cmin_abs))+1);
        conndat_thick_norm(conndat_thick_norm<=0) = 1; % assure no negative or 0 indices
        conndat_thick_norm(~non_nan) = nan;
        conndat_tmap = zeros(numrois); conndat_tmap(triu(true(numrois),1)) = conndat_thick_norm;
        conndat_tmap = conndat_tmap + conndat_tmap';
        % Plot Lines:
        count = 0;
        for ix = 1:numrois-1
            for jxx = ix+1:numrois
                if ~isnan(conndat_cmap(ix,jxx))
                    count = count+1;
                    coords = zeros(2,3);
                    coords(1,:) = mean(roi_patch(ix).p.vertices);
                    coords(2,:) = mean(roi_patch(jxx).p.vertices);
%                     h_line(axes_iter1,count) = plot3(ax1(axes_iter1),coords(:,1),coords(:,2),coords(:,3),'Color',...
%                         rb_cmap(conndat_cmap(ix,jxx),:),'LineWidth',conndat_tmap(ix,jxx),'AlignVertexCenters','on');
                    h_line(axes_iter1,count) = patchline(coords(:,1),coords(:,2),coords(:,3),...
                        'edgecolor',rb_cmap(conndat_cmap(ix,jxx),:),...
                        'LineWidth',conndat_tmap(ix,jxx),'Parent',handles.axes(axes_iter1),'edgealpha',edge_alpha);
                end
            end
        end
    end
end

function change_edges_single_colors(~,~,which_color,initial)
    if nargin>3 && ~isempty(initial)
        % Change Colors:
        for axes_iter1 = 1:n_axes
            for ixx = 1:size(h_line,2)
                set(h_line(axes_iter1,ixx),'edgecolor',initial)
            end
        end
    else
        for ix = 1:8; h_edges_single_colors(ix).Checked = 'off'; end
        h_edges_single_colors(which_color).Checked = 'on';
        if which_color < 8
            % Change Colors:
            for axes_iter1 = 1:n_axes
                for ixx = 1:size(h_line,2)
                    if isgraphics(h_line(axes_iter1,ixx),'patch')
                        set(h_line(axes_iter1,ixx),'EdgeColor',roi_single_colors(which_color,:))
                    end
                end
            end
        else
            if ~ishandle(h_edge_colorbar_fig)
                h_edge_colorbar_fig = figure('menu','none','Name','Choose Color','NumberTitle','off');
                h_edge_colorbar_fig.Position(4) = 70; 
                ax_colorbar = gca; set(ax_colorbar,'Position',[0,0,1,1],'Visible','off','nextplot','add');
                c_map = jet(600); c_image = zeros(100,600,3);
                c_image(:,:,1) = repmat(c_map(:,1)',100,1);
                c_image(:,:,2) = repmat(c_map(:,2)',100,1);
                c_image(:,:,3) = repmat(c_map(:,3)',100,1);
                image(c_image,'Parent',ax_colorbar,'ButtonDownFcn',{@manual_select_edge_color,c_map});
                set(ax_colorbar,'XLim',[0,600],'YLim',[0,100])
            else
                figure(h_edge_colorbar_fig)
            end
        end
    end
end

function manual_select_edge_color(~,~,c_map)
    pt = get(gca,'CurrentPoint');
    xpos = round(pt(1,1));
    use_color = c_map(xpos,:);
    for axes_iter1 = 1:n_axes
        for ix = 1:size(h_line,2)
            set(h_line(axes_iter1,ix),'edgecolor',use_color)
        end
    end
end

function edges_colormap_callback(hObject,~,which_spec)
    switch hObject.Label
        case 'Blue-White-Red'
            rb_cmap = bluewhitered(m,cmin,cmax);
        case 'Red-Blue'
            rb_cmap = redblue(m);
        otherwise
            rb_cmap = eval([hObject.Label,'(',num2str(m),')']);
    end
    for ix = 1:length(edges_color_spec)
        edges_color_spec(ix).Checked = 'off';
    end; edges_color_spec(which_spec).Checked = 'on';
    m = 64;
    change_edges_clim([],[],1)
end

function change_edge_alpha(hObject,~,~)
    edge_alpha = hObject.Value;
    change_edges_clim([],[],1)
end

function keypress_callback(~,which_key,~)
    switch which_key.Key
        case 'uparrow' 
%             rotate_nifti(ax,-5,0,0)
%             nifti_camorbit(h_ax,0,-5)   
            nifti_camorbit(handles.axes,-1,0,'direction','x')  
        case 'downarrow'
%             rotate_nifti(ax,5,0,0)
%             nifti_camorbit(h_ax,0,5)   
            nifti_camorbit(handles.axes,1,0,'direction','x')  
        case 'rightarrow' 
%             rotate_nifti(ax,0,-5,0)
%             nifti_camorbit(h_ax,5,0)   
            nifti_camorbit(handles.axes,-1,0,'direction','z')  
        case 'leftarrow' 
%             rotate_nifti(ax,0,5,0)
%             nifti_camorbit(h_ax,-5,0)  
            nifti_camorbit(handles.axes,1,0,'direction','z')  
        case 'e' 
%             rotate_nifti(h_ax,0,0,-5)
%             roll_nifti(handles.axes,3);
            nifti_camorbit(handles.axes,-1,0,'direction','y')  
        case 'r'
%             rotate_nifti(h_ax,0,0,5)
%             roll_nifti(handles.axes,-3);
            nifti_camorbit(handles.axes,1,0,'direction','y')  
    end
end

end

% nifti_camorbit very slightly modifies the Matlab built-in function
% camorbit.m for better compatibility (Copyright 1984-2012 The MathWorks, Inc.)
function nifti_camorbit(ax, dtheta, dphi, coordsys, direction)
    %CAMORBIT Orbit camera.
    %   CAMORBIT(DTHETA, DPHI)  Orbits (rotates) the camera position
    %   of the current axes around the camera target by the amounts
    %   specified in DTHETA and DPHI (both in degrees). DTHETA is the
    %   horizontal rotation and DPHI is the vertical.
    % 
    %   CAMORBIT(DTHETA, DPHI, coordsys, direction) determines the center
    %   of rotation. Coordsys can be 'data' (the default) or 'camera'.  If
    %   coordsys is 'data' (the default), the camera position rotates
    %   around a line specified by the camera target and direction.
    %   Direction can be 'x', 'y', or 'z' (the default) or [X Y Z].  If
    %   coordsys is 'camera', the rotation is around the camera target
    %   point. 
    %
    %   CAMORBIT(AX, ...) uses axes AX instead of the current axes.
    %
    %   See also CAMDOLLY, CAMPAN, CAMZOOM, CAMROLL.

    %   Copyright 1984-2012 The MathWorks, Inc.

    if nargin>5 || nargin<2
      error(message('MATLAB:camorbit:IncorrectNumberArguments'))
    elseif nargin<5

      if any(ishghandle(ax,'axes'))
        if nargin<3 
          error(message('MATLAB:camorbit:NotEnoughInputs'))
        else
          direction = [0 0 1];
          if nargin==3
        coordsys = 'data';
          end
        end
      else
        if nargin==4
          direction = coordsys;
          coordsys  = dphi;
        elseif nargin==3
          direction = [0 0 1];
          coordsys  = dphi;
        else %nargin==2
          direction = [0 0 1];
          coordsys  = 'data';
        end

        dphi   = dtheta;
        dtheta = ax;
        ax     = gca;
      end

    end
    
%     view_angle = camva(ax); % View Angle is not the culprit for changing zoom
    pos  = get(ax, 'cameraposition' );
    targ = get(ax, 'cameratarget'   );
    dar  = get(ax, 'dataaspectratio');
    up   = get(ax, 'cameraupvector' );

    if ~righthanded(ax), dtheta = -dtheta; end

    [newPos, newUp] = camrotate(pos,targ,dar,up,dtheta,dphi,coordsys,direction);

    if all(isfinite(newPos))  
        set(ax,'CameraPosition',newPos);
    end
    if all(isfinite(newUp))
        set(ax,'CameraUpVector',newUp);
    end
    
%     camva(ax,view_angle)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function val=righthanded(ax)

    dirs=get(ax, {'xdir' 'ydir' 'zdir'}); 
    num=length(find(lower(cat(2,dirs{:}))=='n'));

    val = mod(num,2);
    end

end
