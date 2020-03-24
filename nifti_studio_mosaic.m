function [handles] = nifti_studio_mosaic(varargin)
% NIfTI Studio Mosaic 
%   Plot mosaic of slices of a 3D image (.nii/.nii.gz/.img/.hdr)
%   
% Author:
%   Elliot Layden, The University of Chicago, 2016-19
% 
% Cite: 
%   Layden, E. A. (2019). NIfTI Studio: A Matlab Toolbox for MRI Image 
%   Editing and Visualization.
% 
% Inputs (name-value pair arguments):
% 'background',     The main image; Multiple input types allowed:
%                   1. char/string full-path-to or filename of 
%                   a 3D NIfTI (.nii,.nii.gz) or .img/.hdr pair (specified 
%                   referencing only the .img file)
%                        -note: It is best practice to specify the full 
%                       path to the file, using, 
%                       e.g. "background = fullfile(<path>,<filename>);".
%                   2. struct: may be specified as an already loaded image 
%                   structure as returned by load_nii or load_untouch_nii
%                   from Jimmy Shen's NIfTI Tools
%                   3. matrix: may be specified as a 3D image matrix 
%                   (type: double, single, uint_...)
%                   4. If not specified or is specified as empty [], the 
%                   program will prompt the user to select an image file
% 'overlay',        A second image to be overlaid on 'background (same 
%                   input options as 'background')
% 'slices',         vector specifying which slices to plot (default:
%                   linspace(1,numslices,4))
% 'title',          string denoting title for figure
% 'dimension',      numeric specifying which dimension to extract slices
%                   from (default: 3)
% 'slice_locator',  logical or numeric 1/0, specifying whether or not to
%                   plot 3D slice_locator in last axes (showing depth of
%                   slices)
% 'slice_locator_pos', numeric denoting which axes to plot slice_locator
%                      in; if entered as NaN, will plot the slice_locator 
%                      in a separate figure (default: last axes)
% 'axes_dim',       [1x2 vector] denoting axes dimensions (e.g., [2,3]
%                   would signify plotting the slice mosaic in a figure
%                   with 2 axes row-wise and 3 axes column-wise); note that
%                   product(axe_dim) must == length(slices) if
%                   slice_locator is OFF and == length(slices)+1 if
%                   slice_locator is ON
% 'axes_direction', [numeric: 1 (row-wise) or 2 (column-wise)] denoting which
%                   dimension/direction to proceed in for plottting slices 
%                   within axes (default: 2 (columnwise))
% 'xlim',           specify [vector: lower,upper] xlimits of plots in voxel
%                   units
% 'ylim',           specify [vector: lower,upper] xlimits of plots in voxel
%                   units
% 'zlim',           specify [vector: lower,upper] xlimits of plots in voxel
%                   units
% 'rotate',         [numeric: -180 < rotation < 180] denoting rotation to
%                   apply to all slices in degrees (default: 0) 
% 'slice_labels',   [numeric: 1/0 or logical: true/false, denoting whether
%                   or not to add slice labels
% 'slice_label_phys', 1/0; (0) = use slice #'s, (1) = use physical distance 
%                     from origin
% 'slice_label_pos', [numeric: 1-7]; (1) None, (2) Top Left, (3) Top Middle, 
%                    (4) Top Right, (5) Bottom left, (6) Bottom
%                    Middle, (7) Bottom Right
% 'unit_type',   (1) voxel, (2) physical (default)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Waitbar:
h_wait = waitbar(0, 'Please wait...', 'Name', 'NIfTI Studio Mosaic');
try
    wbch = allchild(h_wait);
    jp = wbch(1).JavaPeer;
    jp.setIndeterminate(1);
catch
end
     
function clearWaitbar
    if exist('h_wait','var') && ~isempty(h_wait) && ishandle(h_wait)
        delete(h_wait)
    end
end

% Identify Function Path and Add Helper Scripts:
script_fullpath = mfilename('fullpath');
[script_path,~,~] = fileparts(script_fullpath);
addpath(genpath(script_path))

roiLimit = 50; % if more ROIs than this specificied, treats as voxel-wise overlay rather than ROIs

% Retrieve name-value pair inputs:
inputs = varargin;
parsed_inputs = struct('background',[],'overlay',[],...
    'slices',[],'title',[],'dimension',3,'slice_locator',0,...
    'slice_locator_pos',[],'axes_dim',[2,5],'axes_direction',2,...
    'customize',0, 'rotate',[],...
    'overlay_alpha',.8,'background_cmap','gray',...
    'background_caxis',[],'overlay_clim',[],'overlay_cmap','jet',...
    'roi_colors',[],'axis_ticks',1,...
    'background_color',zeros(1,3),'mesh_color',[0,0,1],'print',[],'print_res',300,...
    'show_axes',1,'background_thresh',0.00001, 'unit_type', 2, ...
    'physical_units', '', ...
    'xlim',[],'ylim',[],'slice_labels',1,...
    'slice_label_pos',3,'overlay_3D',0, 'origin_slice', 1,'light_axis','x');

poss_input = {'background','overlay','slices','title',...
    'dimension','slice_locator','slice_locator_pos','axes_dim',...
    'axes_direction','customize','rotate','overlay_alpha',...
    'background_cmap','background_caxis','overlay_clim','overlay_cmap',...
    'roi_colors','axis_ticks','background_color','mesh_color','print',...
    'print_res','show_axes','background_thresh','unit_type','physical_units',...
    'xlim','ylim','slice_labels','slice_label_pos',...
    'overlay_3D', 'origin_slice','light_axis'};

m = 200; % # of colormap entries

% Ticks and TickLabels:
xticks = []; yticks = [];
n_ticks_x = 6; n_ticks_y = 10;
h_xticklabels = [];
h_yticklabels = [];
h_units = [0,0];

% Axes tick position
x_tick_position = 2;
y_tick_position = 1;
x_tick_offset = 0.04; % x-tick positioning
y_tick_offset = 0.02; % y-tick positioning

axes_color = ones(1,3);
slice_label_text_color = ones(1,3);

for i = 1:length(poss_input)
    ind = find(strcmp(poss_input{i},inputs));
    if ~isempty(ind)
        input1 = inputs{ind+1};
        parsed_inputs.(poss_input{i}) = input1;
    end
end

unit_type = parsed_inputs.unit_type;
show_axes = parsed_inputs.show_axes;

slice_labels_on = parsed_inputs.slice_labels;
slice_label_pos = parsed_inputs.slice_label_pos;

% Image:
background = parsed_inputs.background;
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
            pixdim = back_img.hdr.dime.pixdim([2,3,4]);
            % Get origin and units
            origin = getOrigin(back_img);
            physical_units = getUnits(back_img);
        case 'struct'
            try
                imdat = background.img;
            catch
                error('If input ''background'' is structure, should contain field ''.img''.')
            end
            pixdim = background.hdr.dime.pixdim([2,3,4]);
            origin = getOrigin(background);
            physical_units = getUnits(background);
        otherwise
            if isnumeric(background)
                imdat = background;
            else
                error('Input ''background'' not recognized.')
            end 
            physical_units = '';
            pixdim = ones(1,3);
    end
else
    load_new_background
end

function physical_units = getUnits(img)
    try
        switch bitand(img.hdr.dime.xyzt_units, 7) % see xform_nii.m, extra_nii_hdr.m
            case 1, physical_units = 'm';
            case 2, physical_units = 'mm';
            case 3, physical_units = 'microns';
            otherwise, physical_units = '';
        end
    catch
        physical_units = '';
        warning('Failed to retrieve voxel units.')
    end
end

function origin = getOrigin(img)
    try
        origin = img.hdr.hist.originator(1:3);
    catch
        origin = img.hdr.dime.dim(2:4)/2;
        warning('Failed to retrieve origin.')
    end
    if any(origin==0)
       origin = img.hdr.dime.dim(2:4)/2;
    end
end

% Load Overlays:
if ~isempty(parsed_inputs.overlay)
    [overlay_dat,~,~] = load_ROI(parsed_inputs.overlay, size(imdat), 'overlay');
    overlay_on = 1;
else
    overlay_dat = [];
    overlay_on = 0;
end

% Permute x- and y-dimensions for display
imdat = permute(imdat,[2,1,3]);
overlay_dat = permute(overlay_dat,[2,1,3]);
pixdim = pixdim([2,1,3]);
origin = origin([2,1,3]);
background_dim = size(imdat);

% Get slice dimensions
dimension = parsed_inputs.dimension; % default 3; never empty
slice_dim = zeros(1,3); slice_pixdim = zeros(1,3); slice_origin = zeros(1,3);
xmax = 0; ymax = 0;
getSliceDimensions(dimension)

set(0,'units','pixels'); 
% screen_res = get(0,'ScreenSize'); % initialize for later resizeFigure

function getSliceDimensions(dim)
    slice_dim = [0,0];
    switch dim
        case 1
            slice_dim = [background_dim(2), background_dim(3), background_dim(1)];
            slice_pixdim = [pixdim(2), pixdim(3), pixdim(1)];
            slice_origin = [origin(2), origin(3), origin(1)];
        case 2
            slice_dim = [background_dim(1), background_dim(3), background_dim(2)];
            slice_pixdim = [pixdim(1), pixdim(3), pixdim(2)];
            slice_origin = [origin(1), origin(3), origin(2)];
        case 3
            slice_dim = [background_dim(2), background_dim(1), background_dim(3)];
            slice_pixdim = [pixdim(2), pixdim(1), pixdim(3)];
            slice_origin = [origin(2), origin(1), origin(3)];
    end
    xmax = slice_dim(1); ymax = slice_dim(2);
end

axes_direction = parsed_inputs.axes_direction; % default 2 (column-wise); never empty
slice_locator_on = parsed_inputs.slice_locator; % default 0 (off); never empty
axes_dim = parsed_inputs.axes_dim; % default [2,5]; never empty (but adjust later if can't accomodate slices)
slices = parsed_inputs.slices; % no default; can be empty
numslices = [];

% Check Slices Specified:
findSlices

% Utility for automatically calculating slices
function findSlices
    % Number of slices (depends on slice_locator)
    if slice_locator_on
        numslices = prod(axes_dim) - 1;
    else
        numslices = prod(axes_dim);
    end
    % Find slices if not specified
    if isempty(slices)
        slices = unique(round(linspace(.2 * slice_dim(dimension), ...
            .8 * slice_dim(dimension), numslices)));
    end
    % Check slices
    if any(slices > slice_dim(dimension)) || any(slices < 1)
        error(['Slice numbers specified are out of bounds. Slices must be between 1 and ',...
            num2str(slice_dim(dimension)),' for dimension ',num2str(dimension),'.'])
    end
end

% Determine Slice Locator Axes Position:
locator_outside = 0;
if slice_locator_on
    if ~isempty(parsed_inputs.slice_locator_pos) && parsed_inputs.slice_locator_pos <= numslices + 1
        slice_locator_pos = parsed_inputs.slice_locator_pos;
    elseif isnan(parsed_inputs.slice_locator_pos)
        locator_outside = 1;
    else
        slice_locator_pos = length(slices)+1; % last axes
    end
    n_axes = numslices+1;
else
    n_axes = numslices;
end

% Adjust axes dimensions if unable to accomodate slices:
if prod(axes_dim) < n_axes
    if slice_locator_on
        warning('The product of input ''axes_dim'' must be greater than or equal to #slices + 1 if the slice_locator is on. Ignoring input...')
    else
        warning('The product of input ''axes_dim'' must be greater than or equal to #slices when the slice_locator is off. Ignoring input...')
    end
    if n_axes<=20 && n_axes>15
        axes_dim = [4,5];
    elseif n_axes<=15 && n_axes>10
        axes_dim = [3,5]; 
    elseif n_axes<=10 && n_axes>8
        axes_dim = [2,5];
    elseif n_axes==8
        axes_dim = [2,4];
    elseif n_axes==6 || n_axes==5
        axes_dim = [2,3];
    elseif n_axes==4
        axes_dim = [2,2];
    else
        axes_dim = [1,n_axes];
    end
end

% Get Physical Units:
if ~isempty(parsed_inputs.physical_units) && ischar(parsed_inputs.physical_units)
    physical_units = parsed_inputs.physical_units;
end

%% Initialize Handles & Options:
overlay_alpha = parsed_inputs.overlay_alpha;
background_cmap = parsed_inputs.background_cmap;
% background_caxis = parsed_inputs.background_caxis;
if isstring(parsed_inputs.overlay_cmap)
    overlay_cmap = eval([parsed_inputs.overlay_cmap,'(',num2str(m),')']);
elseif ismatrix(parsed_inputs.overlay_cmap)
    overlay_cmap = parsed_inputs.overlay_cmap;
else
    error('Input ''overlay_cmap'' should either be a string (e.g., ''jet'') or a nColors x 3 (RGB) matrix')
end
overlay_clim = parsed_inputs.overlay_clim;
if isempty(overlay_clim) && overlay_on
    overlay_clim = [min(overlay_dat(:)), max(overlay_dat(:))];
end 
roi_colors = parsed_inputs.roi_colors;
if ~isempty(roi_colors) 
    if all(size(roi_colors)==[length(unique(overlay_dat(:)))-1,3])
        overlay_cmap = roi_colors; 
    else
        warning('Input ''roi_colors'' should be a matrix of size [numrois,3] RGB values.')
    end
end
    
background_color = parsed_inputs.background_color;
background_thresh = parsed_inputs.background_thresh;
opacity_options = {'Opaque','90%','80%','70%','60%','50%','40%','30%','20%','10%','Invisible'};
imdat(isnan(imdat)) = background_thresh;
if overlay_on; overlay_dat(isnan(overlay_dat)) = background_thresh; end
mesh_color = parsed_inputs.mesh_color;
scroll_zoom_equiv = [1,.9,.8,.7,.6,.5,.4,.3,.2,.1];
scroll_count = 0;
xmin = 1; ymin = 1;
first_click = [];
drag_click = [];
curr_ax = 28.3938;

handles = struct('figure',[],'title',[],'background_axes',[],...
    'background_images',[],'overlay_axes',[],'overlay_images',[],...
    'slice_locator',[],'h_slice_labels',[],'h_slice_labels_rect',[],'mesh',[]);
handles.background_axes = repmat(10.28173,1,n_axes);
handles.overlay_axes = repmat(10.28173,1,n_axes);
handles.background_images = repmat(10.28173,1,n_axes);
handles.overlay_images = repmat(10.28173,1,n_axes);

ax_pos = zeros(n_axes,4);
slice_distances = zeros(1,numslices);

% Initialize slice data containers
background_slice_data = zeros(slice_dim(1),slice_dim(2),numslices);
if overlay_on
    overlay_slice_data = zeros(slice_dim(1),slice_dim(2),numslices);
end

h_colorbar_fig = 28.332737;
h_overlay_colorbar_fig = 37.163224;
h_new_plot_fig = 382.4747;
background_colors = [1,1,1;0,0,0;.2,.6,.7;0,0,1;1,0,0;0,1,0;0,1,1;1,0,1;1,1,0];
roi_single_colors = [1,0,0;0,0,1;0,1,0;0,1,1;1,0,1;1,1,0;0,0,0];
annotation_background_colors = [0,0,0;1,0,0;0,0,1;0,1,0;0,1,1;1,0,1;1,1,0;0,0,0;1,1,1];
first_cmap_change = 1;

% roi_cmap_options = {'jet','hot','cool','hsv','bone','colorcube','copper',...
%             'spring','summer','winter','pink','gray'};

% New Plot Handles:
slices_spec = 35.38848;
dim_bg = 35.38848;
h_dim = repmat(35.38848,1,3);
num_rows_spec = 35.38848;
num_cols_spec = 35.38848;
tile_bg = 35.38848;
h_tile_dir = repmat(35.38848,1,2);
slice_locator_spec = 35.38848;
slice_locator_row_spec = 35.38848;
slice_locator_col_spec = 35.38848;
handles.h_slice_labels = repmat(371.38389,1,numslices);
        
%% FIGURE AND FILE MENUS:
handles.figure = figure('menubar','none','color',background_color,'numbertitle',...
    'off','name','','units','norm','Position',[0,0,1,1],'Pointer','hand');  % [.15,.099,.73,.802]

% % Try to maximize
% if verLessThan('matlab','9.4') % R2018a introduced WindowState feature
%     try 
%         jFrame = get(handles.figure, 'JavaFrame');
%         jFrame.setMaximized(1);   
%     catch
%     end
% else 
%     set(handles.figure, 'WindowState','maximized');
% end

set(handles.figure,'WindowScrollWheelFcn',@scroll_zoom_callback);
set(handles.figure,'WindowButtonDownFcn',@cursor_click_callback);

% FILE MENUS:
color_spectrums = {'Blue-White-Red', 'Red-Blue', 'jet', 'hot', 'cool', ...
    'hsv', 'bone', 'colorcube', 'copper', 'spring', 'summer', ...
    'winter', 'pink', 'gray'};

file_menu = uimenu(handles.figure,'Label','File');
    uimenu(file_menu,'Label','Save Figure','Callback',@save_figure_callback);
    uimenu(file_menu,'Label','Print','Callback',{@print_callback,0});
new_plot_menu = uimenu(handles.figure,'Label','New Plot');
    uimenu(new_plot_menu, 'Label','Specify...','Callback',@new_plot)
tools_menu = uimenu(handles.figure,'Label','Tools');
    uimenu(tools_menu,'Label','Rotate Slices','Callback',@rotate_slices);
    h_tools(1) = uimenu(tools_menu,'Label','Pan','Checked','on','Callback',{@change_tool,1});
    h_tools(2) = uimenu(tools_menu,'Label','Zoom','Checked','on','Callback',{@change_tool,2});
    uimenu(tools_menu,'Label','Revert','Callback',@revert_view);
display_menu = uimenu(handles.figure,'Label','Display');
    measurements_menu = uimenu(display_menu,'Label','Measurements');
            h_units(1) = uimenu(measurements_menu,'Label','Voxels','Checked','off','Callback',{@change_units,1});
            h_units(2) = uimenu(measurements_menu,'Label','Physical','Checked','on','Callback',{@change_units,2});  
    background_options_menu = uimenu(display_menu,'Label','Background');
        background_color_menu = uimenu(background_options_menu, 'Label', 'Color');
            h_background(1) = uimenu(background_color_menu,'Label','White','Callback',{@change_background,1});
            h_background(2) = uimenu(background_color_menu,'Label','Black','Checked','on','Callback',{@change_background,2});
            h_background(3) = uimenu(background_color_menu,'Label','Aqua','Callback',{@change_background,3});
            h_background(4) = uimenu(background_color_menu,'Label','Blue','Callback',{@change_background,4});
            h_background(5) = uimenu(background_color_menu,'Label','Red','Callback',{@change_background,5});
            h_background(6) = uimenu(background_color_menu,'Label','Green','Callback',{@change_background,6});
            h_background(7) = uimenu(background_color_menu,'Label','Cyan','Callback',{@change_background,7});
            h_background(8) = uimenu(background_color_menu,'Label','Magenta','Callback',{@change_background,8});
            h_background(9) = uimenu(background_color_menu,'Label','Yellow','Callback',{@change_background,9});
            h_background(10) = uimenu(background_color_menu,'Label','More...','Callback',{@change_background,10});
        uimenu(background_options_menu,'Label','Intensity Threshold','Callback',@change_background_thresh);
    h_axes_menu = uimenu(display_menu,'Label','Axes');
        h_show_axes = uimenu(h_axes_menu,'Label','Visible', ...
            'Checked','on','Callback', @change_show_axes);
        if ~parsed_inputs.show_axes
            set(h_show_axes, 'Checked','off');
        end
        h_axes_position = uimenu(h_axes_menu, 'Label','Position');
            h_x_axis = uimenu(h_axes_position, 'Label','X');
                h_x_axis_position(1) = uimenu(h_x_axis, 'Label','Top', ...
                    'Checked','on','Callback', {@change_x_position, 1});
                h_x_axis_position(2) = uimenu(h_x_axis, 'Label','Bottom', ...
                    'Checked','off','Callback', {@change_x_position, 2});
            h_y_axis = uimenu(h_axes_position, 'Label','Y');
                h_y_axis_position(1) = uimenu(h_y_axis, 'Label','Left', ...
                    'Checked','on','Callback', {@change_y_position, 1});
                h_y_axis_position(2) = uimenu(h_y_axis, 'Label','Right', ...
                    'Checked','off','Callback', {@change_y_position, 2});
        menu_axis_color = uimenu(h_axes_menu,'Label','Color');
            h_ax_color(1) = uimenu(menu_axis_color,'Label','Grey','Checked','off','Callback',{@change_axes_color,1});
            h_ax_color(2) = uimenu(menu_axis_color,'Label','Black','Checked','off','Callback',{@change_axes_color,2});
            h_ax_color(3) = uimenu(menu_axis_color,'Label','White','Checked','on','Callback',{@change_axes_color,3});
        axes_limits_menu = uimenu(h_axes_menu,'Label','Limits');
            uimenu(axes_limits_menu,'Label','XLim','Callback', @change_x_limits);
            uimenu(axes_limits_menu,'Label','YLim','Callback', @change_y_limits);
    slice_labels_menu = uimenu(display_menu,'Label','Slice Labels');
    h_slice_labels_position_menu = uimenu(slice_labels_menu,'Label','Position');
        h_slice_label_pos_menu(1) = uimenu(h_slice_labels_position_menu,'Label','None','Callback',{@change_slice_label_position,1});
        h_slice_label_pos_menu(2) = uimenu(h_slice_labels_position_menu,'Label','Top Left','Callback',{@change_slice_label_position,2});
        h_slice_label_pos_menu(3) = uimenu(h_slice_labels_position_menu,'Label','Top Middle','Callback',{@change_slice_label_position,3});
        h_slice_label_pos_menu(4) = uimenu(h_slice_labels_position_menu,'Label','Top Right','Callback',{@change_slice_label_position,4});
        h_slice_label_pos_menu(5) = uimenu(h_slice_labels_position_menu,'Label','Bottom Left','Callback',{@change_slice_label_position,5});
        h_slice_label_pos_menu(6) = uimenu(h_slice_labels_position_menu,'Label','Bottom Middle','Callback',{@change_slice_label_position,6});
        h_slice_label_pos_menu(7) = uimenu(h_slice_labels_position_menu,'Label','Bottom Right','Callback',{@change_slice_label_position,7});
        set(h_slice_label_pos_menu(slice_label_pos),'Checked','on')
    slice_labels_color_menu = uimenu(slice_labels_menu,'Label','Color');    
    slice_label_txt_color = uimenu(slice_labels_color_menu,'Label','Text');
        h_slice_label_color(1) = uimenu(slice_label_txt_color,'Label','Grey','Checked','off','Callback',{@change_label_txt_color,1});
        h_slice_label_color(2) = uimenu(slice_label_txt_color,'Label','Black','Checked','on','Callback',{@change_label_txt_color,2});
        h_slice_label_color(3) = uimenu(slice_label_txt_color,'Label','White','Checked','off','Callback',{@change_label_txt_color,3});
    h_slice_labels_background_menu = uimenu(slice_labels_color_menu,'Label','Background');
        h_slice_label_background_color(1) = uimenu(h_slice_labels_background_menu,'Label','None','Checked','on','Callback',{@change_slice_labels_background,1});
        h_slice_label_background_color(2) = uimenu(h_slice_labels_background_menu,'Label','Red','Callback',{@change_slice_labels_background,2});
        h_slice_label_background_color(3) = uimenu(h_slice_labels_background_menu,'Label','Blue','Callback',{@change_slice_labels_background,3});
        h_slice_label_background_color(4) = uimenu(h_slice_labels_background_menu,'Label','Green','Callback',{@change_slice_labels_background,4});
        h_slice_label_background_color(5) = uimenu(h_slice_labels_background_menu,'Label','Cyan','Callback',{@change_slice_labels_background,5});
        h_slice_label_background_color(6) = uimenu(h_slice_labels_background_menu,'Label','Magenta','Callback',{@change_slice_labels_background,6});
        h_slice_label_background_color(7) = uimenu(h_slice_labels_background_menu,'Label','Yellow','Callback',{@change_slice_labels_background,7});
        h_slice_label_background_color(8) = uimenu(h_slice_labels_background_menu,'Label','Black','Callback',{@change_slice_labels_background,8});
        h_slice_label_background_color(9) = uimenu(h_slice_labels_background_menu,'Label','White','Callback',{@change_slice_labels_background,9});
        h_slice_label_background_color(10) = uimenu(h_slice_labels_background_menu,'Label','More...','Callback',{@change_slice_labels_background,10});
% Add Overlay Menu if On:
if overlay_on
    overlay_menu = uimenu(handles.figure,'Label','Overlay');
        uimenu(overlay_menu,'Label','Color Limits','Callback',@change_overlay_clim);
        overlay_alpha_menu = uimenu(overlay_menu,'Label','Opacity');
            overlay_alpha_submenus = zeros(1,11);
            for i = 1:11
                overlay_alpha_submenus(i) = uimenu(overlay_alpha_menu,...
                    'Label',opacity_options{i},'Callback',{@overlay_opacity_callback,i});
            end
        overlay_colors_menu = uimenu(overlay_menu,'Label','Color Scheme');
            overlay_color_spectrum_menu = uimenu(overlay_colors_menu,'Label','Color Spectrum');
                h_overlay_color_spec = zeros(1,length(color_spectrums));
                for i = 1:length(color_spectrums)
                    h_overlay_color_spec(i) = uimenu(overlay_color_spectrum_menu,...
                        'Label',color_spectrums{i}, ...
                        'Callback',{@overlay_colormap_callback,i});
                end
            overlay_single_colors_menu = uimenu(overlay_colors_menu,'Label','Single Color');
                h_overlay_single_colors(1) = uimenu(overlay_single_colors_menu,'Label','Red','Callback',{@change_overlay_single_colors,1});
                h_overlay_single_colors(2) = uimenu(overlay_single_colors_menu,'Label','Blue','Callback',{@change_overlay_single_colors,2});
                h_overlay_single_colors(3) = uimenu(overlay_single_colors_menu,'Label','Green','Callback',{@change_overlay_single_colors,3});
                h_overlay_single_colors(4) = uimenu(overlay_single_colors_menu,'Label','Cyan','Callback',{@change_overlay_single_colors,4});
                h_overlay_single_colors(5) = uimenu(overlay_single_colors_menu,'Label','Magenta','Callback',{@change_overlay_single_colors,5});
                h_overlay_single_colors(6) = uimenu(overlay_single_colors_menu,'Label','Yellow','Callback',{@change_overlay_single_colors,6});
                h_overlay_single_colors(7) = uimenu(overlay_single_colors_menu,'Label','Black','Callback',{@change_overlay_single_colors,7});
                h_overlay_single_colors(8) = uimenu(overlay_single_colors_menu,'Label','More...','Callback',{@change_overlay_single_colors,8});
end
% Add Slice Locator Menu:
if slice_locator_on
    slice_locator_menu = uimenu(handles.figure,'Label','Slice Locator');
        mesh_options = uimenu(slice_locator_menu,'Label','Slice Mesh');
        mesh_face_alpha_menu = uimenu(mesh_options,'Label','Face Opacity');
            mesh_face_alpha_submenus = zeros(1,11);
            for i = 1:11
                mesh_face_alpha_submenus(i) = uimenu(mesh_face_alpha_menu,...
                    'Label',opacity_options{i},'Callback',{@mesh_face_opacity_callback,i});
            end
            set(mesh_face_alpha_submenus(8),'Checked','on')
        mesh_edge_alpha_menu = uimenu(mesh_options,'Label','Edge Opacity');
            mesh_edge_alpha_submenus = zeros(1,11);
            for i = 1:11
                mesh_edge_alpha_submenus(i) = uimenu(mesh_edge_alpha_menu,'Label',opacity_options{i},'Callback',{@mesh_edge_opacity_callback,i});
            end
            set(mesh_edge_alpha_submenus(11),'Checked','on')
end
% End Menus

% Allow advanced user customization
if parsed_inputs.customize
    new_plot([], [], 1)
else
    plot_mosaic
end

% Adjust colormap:
if overlay_on
    if any(strcmp(parsed_inputs.overlay_cmap, color_spectrums))
        overlay_colormap_callback(parsed_inputs.overlay_cmap, [], ...
            find(strcmp(parsed_inputs.overlay_cmap, color_spectrums)))
    else
        warning('Input ''overlay_cmap'' not found as option.')
    end
end

% Change unit of measure if specified:
change_units([],[], unit_type)

% Adjust axes position
change_x_position([],[], 2) % default: bottom
change_y_position([],[], 1) % default: left

% Delete waitbar
clearWaitbar

%% Begin plotting
function plot_mosaic
    % Determine Axes Positions
    row_s = 1/axes_dim(1); 
    col_s = 1/axes_dim(2);
    row_l = 1:-row_s:0; row_l = row_l(2:end);
    col_l = 0:col_s:1; col_l = col_l(1:end-1);
    count = 0;
    if axes_direction==2
        for row = 1:axes_dim(1)
            for col = 1:axes_dim(2)
                if count<numslices
                    count = count + 1;
                    ax_pos(count,:) = [col_l(col),row_l(row),col_s,row_s];
                end
            end
        end
    elseif axes_direction==1
        for col = 1:axes_dim(2)
            for row = 1:axes_dim(1)
                if count<numslices
                    count = count + 1;
                    ax_pos(count,:) = [col_l(col),row_l(row),col_s,row_s];
                end
            end
        end
    end
    % ADD SLICE LOCATOR
    if slice_locator_on
        ax_pos(end,:) = [col_l(end),row_l(end),col_s,row_s];
        handles.slice_locator = axes('parent',handles.figure,'Box','off',...
            'position',ax_pos(end,:),'NextPlot','add','Visible','off',...
            'XTick',[],'YTick',[]);
        slice_pos_dat = zeros(size(imdat)); 
        for ix = 1:numslices
            switch dimension
                case 1, slice_pos_dat(slices(ix),:,:) = ix;
                case 2, slice_pos_dat(:,slices(ix),:) = ix;
                case 3, slice_pos_dat(:,:,slices(ix)) = ix;
            end
        end
        % Add barrier between Mosaic menus and nifti_studio_3D menus:
        uimenu(handles.figure,'Label',['                                                   ',...
    '                                                                                   ',...
    '                                   ']);
        if parsed_inputs.overlay_3D 
            if locator_outside
                if numel(unique(overlay_dat(:))) > roiLimit % don't depict ROIs
                    [~] = nifti_studio_3D('background',imdat,...
                        'show_axes',0,'view_spec',[125,-5],'view_angle',7.5,...
                        'background_alpha',.36,'show_grid',0,'menu_on',1,...
                        'xlim',parsed_inputs.ylim,'ylim',parsed_inputs.xlim,...
                        'light_axis', parsed_inputs.light_axis, 'background_color',zeros(1,3));
                else % depict ROIs
                    [~] = nifti_studio_3D('background',imdat,'ROI',overlay_dat,...
                        'roi_colors',roi_colors,'show_axes',0,...
                        'view_spec',[125,-5],'view_angle',7.5,'background_alpha',.36,...
                        'roi_alpha',1,'roi_type','cluster','show_grid',0,'menu_on',1,...
                        'xlim',parsed_inputs.ylim,'ylim',parsed_inputs.xlim, ...
                        'light_axis', parsed_inputs.light_axis, 'background_color',zeros(1,3));
                end
            else
                if numel(unique(overlay_dat(:))) > roiLimit % don't depict ROIs
                    [~] = nifti_studio_3D('background',imdat,'insert_axes',handles.slice_locator,...
                        'show_axes',0,'view_spec',[125,-5],'view_angle',7.5,...
                        'background_alpha',.36,'show_grid',0,'menu_on',1,...
                        'xlim',parsed_inputs.ylim,'ylim',parsed_inputs.xlim,...
                        'light_axis', parsed_inputs.light_axis, 'background_color',zeros(1,3));
                else
                    [~] = nifti_studio_3D('background',imdat,'ROI',overlay_dat,...
                        'roi_colors',roi_colors,'insert_axes',handles.slice_locator,'show_axes',0,...
                        'view_spec',[125,-5],'view_angle',7.5,'background_alpha',.36,...
                        'roi_alpha',1,'roi_type','cluster','show_grid',0,'menu_on',1,...
                        'xlim',parsed_inputs.ylim,'ylim',parsed_inputs.xlim,...
                        'light_axis', parsed_inputs.light_axis, 'background_color',zeros(1,3));
                end
            end
        else
            if locator_outside
                [~] = nifti_studio_3D('background',imdat,...
                    'show_axes',0,...
                    'background_alpha',1,'roi_type','cluster','show_grid',0,'menu_on',1,...
                    'xlim',parsed_inputs.ylim,'ylim',parsed_inputs.xlim,...
                    'light_axis', parsed_inputs.light_axis, 'background_color',zeros(1,3));
                    %'view_spec',[125,-5],'cam_roll',-75,'view_angle',6.4,...
            else
                [~] = nifti_studio_3D('background',imdat,...
                    'insert_axes',handles.slice_locator,'show_axes',0,...
                    'background_alpha',1,'roi_type','cluster','show_grid',0,'menu_on',1,...
                    'xlim',parsed_inputs.ylim,'ylim',parsed_inputs.xlim,...
                    'light_axis', parsed_inputs.light_axis, 'background_color',zeros(1,3));
                    %... % 'view_spec',[125,-5],'cam_roll',-75,'view_angle',6.4,

            end
        end
%         rotate3d off
        % Add mesh for slice locations:
        mesh_step = 15;
        mesh_coords_x = 1:mesh_step:slice_dim(1); 
        mesh_coords_y = 1:mesh_step:slice_dim(2);
        mesh_CData = zeros(length(mesh_coords_y),length(mesh_coords_x),3);
        for ix = 1:length(mesh_coords_y)
            for iy = 1:length(mesh_coords_x)
                mesh_CData(ix,iy,:) = mesh_color;
            end
        end
        handles.mesh = zeros(1,numslices);
    end
    
    % Create Axes and Plot Mosaic:
    for ix = 1:numslices
        handles.background_axes(ix) = axes('parent',handles.figure,...
            'position',ax_pos(ix,:),'NextPlot','add','units','normalized',...
            'color', background_color);
        colormap(handles.background_axes(ix),background_cmap);
        if ~parsed_inputs.show_axes
            set(handles.background_axes(ix),'Visible','off');
        end
        switch dimension
            case 1, background_slice = squeeze(imdat(slices(ix),:,:))';
                set(handles.background_axes(ix),'Xdir','normal')
            case 2, background_slice = squeeze(imdat(:,slices(ix),:))';
                set(handles.background_axes(ix),'Xdir','reverse')
            case 3, background_slice = squeeze(imdat(:,:,slices(ix)));
                set(handles.background_axes(ix),'Xdir','normal')
        end
        set(handles.background_axes(ix),...
            'XLim',[1, slice_dim(1)],...
            'YLim',[1, slice_dim(2)]);
        if all(size(background_slice_data(:,:,ix)) == size(background_slice))
            background_slice_data(:,:,ix) = background_slice;
        else 
            background_slice_data(:,:,ix) = background_slice';
        end
        alpha_background = zeros(size(background_slice));
        alpha_background(background_slice>=background_thresh) = 1; 
        handles.background_images(ix) = imagesc(background_slice,'Parent',...
            handles.background_axes(ix),'AlphaDataMapping','none',...
            'AlphaData',alpha_background);
        if overlay_on
            handles.overlay_axes(ix) = axes('parent',handles.figure,...
                'position',ax_pos(ix,:),'NextPlot','add',...
                'Visible','off','units','normalized');
            colormap(handles.overlay_axes(ix),overlay_cmap);
            switch dimension
                case 1
                    overlay_slice = squeeze(overlay_dat(slices(ix),:,:))';
                    set(handles.overlay_axes(ix),'Xdir','normal')
                case 2
                    overlay_slice = squeeze(overlay_dat(:,slices(ix),:))';
                    set(handles.overlay_axes(ix),'Xdir','reverse')
                case 3
                    overlay_slice = squeeze(overlay_dat(:,:,slices(ix)));
                    set(handles.overlay_axes(ix),'Xdir','normal')
            end
            set(handles.overlay_axes(ix),...
                'XLim',[1, slice_dim(1)],...
                'YLim',[1, slice_dim(2)]);
            alpha_slice = zeros(size(overlay_slice));
            alpha_slice(overlay_slice>overlay_clim(1)) = overlay_alpha;
            alpha_slice(overlay_slice>overlay_clim(2)) = 0;
            alpha_slice(overlay_slice == 0) = 0; % set zeros to transparent
%             alpha_slice(overlay_slice>overlay_clim(2)) = 0;
            
            if all(size(overlay_slice_data(:,:,ix)) == size(overlay_slice))
                overlay_slice_data(:,:,ix) = overlay_slice;
            else 
                overlay_slice_data(:,:,ix) = overlay_slice';
            end
        
            % Scale Slice Data:
            if ~isempty(roi_colors)
                handles.overlay_images(ix) = image(overlay_slice,'Parent',...
                    handles.overlay_axes(ix),'AlphaDataMapping','none',...
                    'AlphaData',alpha_slice,'CDataMapping','direct');
            else
                m = 100; 
                overlay_vec = overlay_slice(:);
%                 ind_nan = ((overlay_vec<=overlay_clim(1))+(overlay_vec>overlay_clim(2)))>0;
                ind_nan = overlay_vec<=overlay_clim(1);
                overlay_vec = min(m,round((m-1).*(overlay_vec-overlay_clim(1))./(overlay_clim(2)-overlay_clim(1)))+1);
                overlay_vec(overlay_vec<=0) = 1; % assure no negative or 0 indices
                overlay_vec(ind_nan) = nan; 
                overlay_slice(:) = overlay_vec;
                handles.overlay_images(ix) = image(overlay_slice,'Parent',...
                    handles.overlay_axes(ix),'AlphaDataMapping','none',...
                    'AlphaData',alpha_slice,'CDataMapping','direct');
            end
        end
        if slice_locator_on
            [mesh_X,mesh_Y,mesh_Z] = meshgrid(mesh_coords_x,mesh_coords_y,slices(ix));
            handles.mesh(ix) = mesh(handles.slice_locator,mesh_X,mesh_Y,mesh_Z,'CData',mesh_CData,...
                'FaceLighting','gouraud','AlphaData',.3,'AlphaDataMapping',...
                'direct','FaceAlpha',.3,'EdgeAlpha',0,'FaceColor','Interp');
        end
    end
    
    xmin = 1; xmax = slice_dim(1);
    ymin = 1; ymax = slice_dim(2);
    
    % Toggle show_axes
    change_show_axes([], [])
    
    % Update slice labels
    change_slice_label_position([], [], slice_label_pos)   
    
    % Update axes color:
    if all(background_color==0)
        change_axes_color([], [], 3) % white
        change_label_txt_color([], [], 3)
    else
        change_axes_color([], [], 1) % grey
        change_label_txt_color([], [], 1)
    end
    
    % Resize figure based on ratio of x- and y- axes
    resizeFigure
    
    % TickLabels:
    adjust_tick_labels
end

rotate3d off
if overlay_on && ishandle(handles.overlay_axes(1))
    axes(handles.overlay_axes(1));
elseif ishandle(handles.background_axes(1))
    axes(handles.background_axes(1));
end

% Rotate if Requested:
if ~isempty(parsed_inputs.rotate) && isnumeric(parsed_inputs.rotate)
    rotate_slices([], [], parsed_inputs.rotate);
end

% Change X-Lim and/or Y-Lim if requested:
if ~isempty(parsed_inputs.xlim)
    change_x_limits([],[], parsed_inputs.xlim);
end
if ~isempty(parsed_inputs.ylim)
    change_y_limits([],[], parsed_inputs.ylim);
end

% Add Slice Labels if Requested:
if ~slice_labels_on
    slice_label_pos = 1;
end
change_slice_label_position([], [], slice_label_pos)

% Print from command line:
if ~isempty(parsed_inputs.print) && ischar(parsed_inputs.print)
    print_callback([],[],1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE ACTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scroll_zoom_callback(~, eventdata, ~)
    % Get mouse location
    curr_ax = neuroimage_overobj('axes'); % find which axes mouse is over
    if isempty(curr_ax); return; end
    pt = get(curr_ax,'CurrentPoint'); 
    if isempty(pt); return; end % returns [] if outside axes
    pt = round(pt(1,1:2));
    % Check if Image Click
    if pt(2) <= slice_dim(1) && pt(2) > 0 && pt(1) <= slice_dim(2) && pt(1) > 0
        scroll_count = min(max(scroll_count + -eventdata.VerticalScrollCount,0),9);
        if scroll_count ~= 0
            zoom_factor = scroll_zoom_equiv(scroll_count+1);
            xlength = zoom_factor * slice_dim(1); 
            ylength = zoom_factor * slice_dim(2);
            xmin = max(round(pt(1)-.5*xlength),1);
            xmax = min(round(pt(1)+.5*xlength), slice_dim(1));
            ymin = max(round(pt(2)-.5*ylength),1);
            ymax = min(round(pt(2)+.5*ylength), slice_dim(2));
        else
            xmin = 1; xmax = slice_dim(1);
            ymin = 1; ymax = slice_dim(2);
        end
        for ix = 1:numslices
            set(handles.background_axes(ix),'XLim',[xmin,xmax],'YLim',[ymin,ymax]);
            if overlay_on && isgraphics(handles.overlay_axes(ix),'axes')
                set(handles.overlay_axes(ix),'XLim',[xmin,xmax],'YLim',[ymin,ymax]);
            end
        end
    end
    adjust_tick_labels
end

function cursor_click_callback(~,~,~)
    curr_ax = neuroimage_overobj('axes'); % find which axes mouse is over
    if isempty(curr_ax); return; end
    pt = get(curr_ax,'CurrentPoint'); 
    if isempty(pt); return; end % returns [] if outside axes
    first_click = round(pt(1,1:2));
    % Set new callbacks:
    set(handles.figure,'WindowButtonMotionFcn', @cursor_motion_callback);
    set(handles.figure,'WindowButtonUpFcn', @cursor_unclick_callback);
end

function cursor_motion_callback(~,~,~)
    % Get mouse location
    curr_ax = neuroimage_overobj('axes'); % find which axes mouse is over
    if isempty(curr_ax); return; end
    pt = get(curr_ax,'CurrentPoint'); 
    if isempty(pt); return; end % returns [] if outside axes
    drag_click = round(pt(1,1:2));
    pan_dir = drag_click-first_click;
    pan_dir = -1.*pan_dir;
    xmin = xmin + pan_dir(1); 
    xmax = xmax + pan_dir(1);
    ymin = ymin + pan_dir(2); 
    ymax = ymax + pan_dir(2);
    for ix = 1:numslices
        set(handles.background_axes(ix),'XLim',[xmin,xmax],'YLim',[ymin,ymax]);
        if overlay_on && isgraphics(handles.overlay_axes(ix),'axes')
            set(handles.overlay_axes(ix),'XLim',[xmin,xmax],'YLim',[ymin,ymax]);
        end
    end
    if overlay_on
        axes(handles.overlay_axes(1))
    end
    adjust_tick_labels
end

function cursor_unclick_callback(~,~,~)
    set(handles.figure,'WindowButtonMotionFcn',[]);
    set(handles.figure,'WindowButtonUpFcn',[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILE MENU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save Figure Function:
function save_figure_callback(~,~,~)
    % Identify File Extension:
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW PLOT MENU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function new_plot(~,~,~)
    if ~ishandle(h_new_plot_fig) 
        if ~isempty(imdat) && nargin < 3
            load_q = questdlg('Load new image?','Load New');
            if strcmp(load_q,'Yes')
                status = load_new_background;
                if ~status % was cancelled
                    return
                end
            elseif strcmp(load_q,'No')
            else
                return
            end
        end
        font_color = [.873,.546,.347];
        h_new_plot_fig = figure('menubar','none','color',zeros(1,3),'numbertitle',...
            'off','name','New Plot Settings','units','norm','Position',[.3,.3,.4,.4]);  % .25,.16,.51,.69
        slices_txt = uicontrol('Style','text','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.03,.82,.42,.14],'FontName',...
            'Helvetica','FontSize',12,'BackgroundColor',zeros(1,3),'String',...
            sprintf('Enter slice numbers \nseparated by commas:'),...
            'ForegroundColor',font_color,'FontWeight','bold','HorizontalAlignment','left'); %#ok
        slices_spec = uicontrol('Style','edit','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.51,.838,.46,.1],'FontName',...
            'Helvetica','FontSize',12);
        slice_dim_txt = uicontrol('Style','text','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.03,.69,.3,.07],'FontName',...
            'Helvetica','FontSize',12,'BackgroundColor',zeros(1,3),'String',...
            'Slice Dimension:','ForegroundColor',font_color,'FontWeight',...
            'bold','HorizontalAlignment','left'); %#ok
        xpos = linspace(.07,.85,3);
        % Slice dimension button group
        dim_bg = uibuttongroup('Visible','on','Units','normalized','Position',...
            [.4,.678,.57,.1],'BackgroundColor',zeros(1,3),'BorderType','etchedin',... 
            'Parent',h_new_plot_fig);
        h_dim(1) = uicontrol(dim_bg,'Style','radiobutton','Units',...
            'normalized','Position',[xpos(1),.23,.12,.6],'HorizontalAlignment','center',...
            'BackgroundColor',zeros(1,3),'String','X','FontName','Helvetica',...
            'FontSize',10,'ForegroundColor',font_color,'FontWeight','bold');
        h_dim(2) = uicontrol(dim_bg,'Style','radiobutton','Units',...
            'normalized','Position',[xpos(2),.23,.12,.6],'HorizontalAlignment','center',...
            'BackgroundColor',zeros(1,3),'String','Y','FontName','Helvetica',...
            'FontSize',10,'ForegroundColor',font_color,'FontWeight','bold');
        h_dim(3) = uicontrol(dim_bg,'Style','radiobutton','Units',...
            'normalized','Position',[xpos(3),.23,.12,.6],'HorizontalAlignment','center',...
            'BackgroundColor',zeros(1,3),'String','Z','FontName','Helvetica',...
            'FontSize',10,'ForegroundColor',font_color,'FontWeight','bold');
        set(dim_bg,'SelectedObject',h_dim(3))
        mosaic_dim_txt = uicontrol('Style','text','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.03,.535,.35,.07],'FontName',...
            'Helvetica','FontSize',12,'BackgroundColor',zeros(1,3),'String',...
            'Mosaic Dimension:','ForegroundColor',font_color,'FontWeight',...
            'bold','HorizontalAlignment','left'); %#ok
        % Mosaic dimension text boxes
        num_rows_txt = uicontrol('Style','text','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.51,.53,.13,.06],'FontName',...
            'Helvetica','FontSize',10,'BackgroundColor',zeros(1,3),'String',...
            '# Rows:','ForegroundColor',font_color,'FontWeight','bold'); %#ok
        num_cols_txt = uicontrol('Style','text','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.73,.53,.17,.06],'FontName',...
            'Helvetica','FontSize',10,'BackgroundColor',zeros(1,3),'String',...
            '# Columns:','ForegroundColor',font_color,'FontWeight','bold'); %#ok
        % Mosaic dimension edit texts
        num_rows_spec = uicontrol('Style','edit','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.65,.54,.06,.06],'FontName',...
            'Helvetica','FontSize',10);
        num_cols_spec = uicontrol('Style','edit','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.91,.54,.06,.06],'FontName',...
            'Helvetica','FontSize',10);
        % Tile direction (column-wise, row-wise)
        tile_direction_txt = uicontrol('Style','text','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.03,.37,.25,.07],'FontName',...
            'Helvetica','FontSize',12,'BackgroundColor',zeros(1,3),'String',...
            'Tile Direction:','ForegroundColor',font_color,'FontWeight',...
            'bold','HorizontalAlignment','left'); %#ok
        % Tile direction button group
        tile_bg = uibuttongroup('Visible','on','Units','normalized','Position',...
            [.4,.355,.57,.1],'BackgroundColor',zeros(1,3),'BorderType','etchedin',...
            'Parent',h_new_plot_fig);
        h_tile_dir(1) = uicontrol(tile_bg,'Style','radiobutton','Units',...
            'normalized','Position',[.15,.23,.45,.6],'HorizontalAlignment','center',...
            'BackgroundColor',zeros(1,3),'String','Column-wise','FontName',...
            'Helvetica','FontSize',10,'ForegroundColor',font_color,'FontWeight','bold');
        h_tile_dir(2) = uicontrol(tile_bg,'Style','radiobutton','Units',...
            'normalized','Position',[.6,.23,.36,.6],'HorizontalAlignment','center',...
            'BackgroundColor',zeros(1,3),'String','Row-wise','FontName',...
            'Helvetica','FontSize',10,'ForegroundColor',font_color,'FontWeight','bold');
        set(tile_bg,'SelectedObject',h_tile_dir(1))
        slice_locator_txt = uicontrol('Style','text','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.03,.22,.25,.07],'FontName',...
            'Helvetica','FontSize',12,'BackgroundColor',zeros(1,3),'String',...
            'Slice Locator:','ForegroundColor',font_color,'FontWeight',...
            'bold','HorizontalAlignment','left'); %#ok
        slice_locator_spec = uicontrol('Style','Checkbox','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.47,.22,.035,.06],'FontName',...
            'Helvetica','FontSize',10,'BackgroundColor',zeros(1,3),'Value',0);
        % Slice locator position spec text
        slice_locator_row_txt = uicontrol('Style','text','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.54,.21,.1,.06],'FontName',...
            'Helvetica','FontSize',10,'BackgroundColor',zeros(1,3),'String',...
            'Row:','ForegroundColor',font_color,'FontWeight','bold'); %#ok
        slice_locator_col_txt = uicontrol('Style','text','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.75,.21,.15,.06],'FontName',...
            'Helvetica','FontSize',10,'BackgroundColor',zeros(1,3),'String',...
            'Column:','ForegroundColor',font_color,'FontWeight','bold'); %#ok
        % Slice locator spec edit text
        slice_locator_row_spec = uicontrol('Style','edit','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.65,.22,.06,.06],'FontName',...
            'Helvetica','FontSize',10,'BackgroundColor',ones(1,3),'String',...
            num_rows_spec.String);
        slice_locator_col_spec = uicontrol('Style','edit','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.91,.22,.06,.06],'FontName',...
            'Helvetica','FontSize',10,'BackgroundColor',ones(1,3),'String',...
            num_cols_spec.String); 
        % Plot UIbutton
        plot_button = uicontrol('Style','pushbutton','Parent',h_new_plot_fig,...
            'Units','normalized','Position',[.4,.06,.2,.1],'FontName',...
            'Helvetica','FontSize',12,'FontWeight','bold','BackgroundColor',font_color,...
            'String','Plot','ForegroundColor',zeros(1,3),'Callback',@run_new_plot); %#ok
    else
       figure(h_new_plot_fig)
    end
end

function run_new_plot(~,~,~)
    for ix = 1:3
        if get(h_dim(ix),'Value')
            dimension = ix;
            break; 
        end
    end
    
    % Get new slice dimensions
    getSliceDimensions(dimension)
    
    if isempty(slices_spec.String)
        findSlices
    else
        slices = str2double(strsplit(slices_spec.String,','));
    end
    numslices = length(slices);
        
    if ~isempty(num_rows_spec.String)
        axes_dim(1) = str2double(num_rows_spec.String);
    else
        axes_dim(1) = 2;
    end
    if ~isempty(num_cols_spec.String)
        axes_dim(2) = str2double(num_cols_spec.String);
    else 
        axes_dim(2) = 5;
    end
    
    % Slice tile direction toggle
    if get(h_tile_dir(1),'Value')
        axes_direction = 2;
    elseif get(h_tile_dir(2),'Value')
        axes_direction = 1;
    end
    
    % Slice locator
    slice_locator_on = get(slice_locator_spec,'Value');
    if isempty(slice_locator_row_spec.String) ...
        || isempty(slice_locator_col_spec.String) ...
            || ~isnumeric(str2double(slice_locator_row_spec.String)) ...
                || ~isnumeric(str2double(slice_locator_col_spec.String))
        slice_locator_pos = prod(axes_dim); 
    else
        mat = false(axes_dim); 
        mat(str2double(slice_locator_row_spec.String),str2double(slice_locator_col_spec.String)) = true;
        count = 0;
        if axes_direction==1
            for col = 1:axes_dim(2)
                for row = 1:axes_dim(1)
                    count = count+1;
                    if mat(row,col); break; end
                end
                if mat(row,col); break; end
            end
        elseif axes_direction==2
            for row = 1:axes_dim(1)
                for col = 1:axes_dim(2)
                    count = count+1;
                    if mat(row,col); break; end
                end
                if mat(row,col); break; end
            end
        end
        slice_locator_pos = count;
    end
    
    % Delete old axes ticks:
    if show_axes
        change_show_axes([], [])
    end
    
    % Delete old plot
    for ix = 1:numslices
        if isgraphics(handles.background_axes(ix),'axes')
            delete(handles.background_axes(ix))
        end
        if overlay_on && isgraphics(handles.overlay_axes(ix),'axes')
            delete(handles.overlay_axes(ix))  
        end
    end
    
    if isgraphics(handles.slice_locator,'axes')
        delete(handles.slice_locator)
    end

    % Re-initialize slice data containers
    background_slice_data = zeros(slice_dim(1),slice_dim(2),numslices);
    if overlay_on
        overlay_slice_data = zeros(slice_dim(1),slice_dim(2),numslices);
    end

    plot_mosaic
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOOLS MENU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rotation:
function rotate_slices(~, ~, rotate_degrees)
    if nargin < 3 || isempty(rotate_degrees)
        prompt = {sprintf('Specify rotation degrees:')};
        dlg_title = 'Rotate'; num_lines = [1,35;]; 
        answer = inputdlg(prompt,dlg_title,num_lines);
        if isempty(answer); return; end
        rotate_degrees = str2double(answer{1});
        for ix = 1:numslices
            A = get(handles.background_images(ix),'CData');
            B = imrotate(A,rotate_degrees,'nearest','crop'); % 'nearest','bilinear','bicubic'
            C = get(handles.background_images(ix),'AlphaData');
            D = imrotate(C,rotate_degrees,'nearest','crop'); % 'nearest','bilinear','bicubic'
            set(handles.background_images(ix),'CData',B,'AlphaData',D);
            if overlay_on && isgraphics(handles.overlay_axes(ix),'axes')
                A = get(handles.overlay_images(ix),'CData');
                B = imrotate(A,rotate_degrees,'nearest','crop'); %  'nearest','bilinear','bicubic'
                C = get(handles.overlay_images(ix),'AlphaData');
                D = imrotate(C,rotate_degrees,'nearest','crop'); % 'nearest','bilinear','bicubic'
                set(handles.overlay_images(ix),'CData',B,'AlphaData',D);
            end
        end
    else
        for ix = 1:numslices
                A = get(handles.background_images(ix),'CData');
                B = imrotate(A,rotate_degrees,'nearest','crop'); % 'nearest','bilinear','bicubic'
                C = get(handles.background_images(ix),'AlphaData');
                D = imrotate(C,rotate_degrees,'nearest','crop'); % 'nearest','bilinear','bicubic'
                set(handles.background_images(ix),'CData',B,'AlphaData',D);
                if overlay_on && isgraphics(handles.overlay_axes(ix),'axes')
                    A = get(handles.overlay_images(ix),'CData');
                    B = imrotate(A,rotate_degrees,'nearest','crop'); %  'nearest','bilinear','bicubic'
                    C = get(handles.overlay_images(ix),'AlphaData');
                    D = imrotate(C,rotate_degrees,'nearest','crop'); % 'nearest','bilinear','bicubic'
                    set(handles.overlay_images(ix),'CData',B,'AlphaData',D);
                end
        end
    end
end

% Pan & Zoom:
function change_tool(~, ~, which_tool)
    axes(handles.background_axes(1))
    rotate3d off; zoom off; pan off;
    switch which_tool
        case 1 % Pan
            if strcmp(h_tools(1).Checked,'on')
                h_tools(1).Checked = 'off';
                set(handles.figure,'WindowButtonDownFcn',[],'Pointer','arrow');
            elseif strcmp(h_tools(1).Checked,'off') 
                h_tools(1).Checked = 'on';
                set(handles.figure,'WindowButtonDownFcn',@cursor_click_callback,'Pointer','hand');
            end
        case 2 % Zoom
            if strcmp(h_tools(2).Checked,'on')
                h_tools(2).Checked = 'off';
                set(handles.figure,'WindowScrollWheelFcn',[]);
            elseif strcmp(h_tools(2).Checked,'off') 
                h_tools(2).Checked = 'on';
                set(handles.figure,'WindowScrollWheelFcn',@scroll_zoom_callback);
            end
    end
end

function revert_view(~,~,~)
    xmin = 1; xmax = slice_dim(1);
    ymin = 1; ymax = slice_dim(2);
    for ix = 1:numslices
        switch dimension
            case 1
                background_slice = squeeze(imdat(slices(ix),:,:))';
                set(handles.background_axes(ix),'Xdir','normal');
            case 2
                background_slice = squeeze(imdat(:,slices(ix),:))';
                set(handles.background_axes(ix),'Xdir','reverse');
            case 3
                background_slice = squeeze(imdat(:,:,slices(ix)));
                set(handles.background_axes(ix),'Xdir','normal');
        end
        set(handles.background_axes(ix),...
            'XLim',[xmin, slice_dim(1)],...
            'YLim',[ymin, slice_dim(2)]);
        
        alpha_background = zeros(size(background_slice));
        alpha_background(background_slice>background_thresh) = 1; 
        set(handles.background_images(ix),'CData',background_slice,...
            'AlphaDataMapping','none','AlphaData',alpha_background);
        
        if overlay_on
            switch dimension
                case 1
                    overlay_slice = squeeze(overlay_dat(slices(ix),:,:))';
                    set(handles.overlay_axes(ix),'Xdir','normal')
                case 2
                    overlay_slice = squeeze(overlay_dat(:,slices(ix),:))';
                    set(handles.overlay_axes(ix),'Xdir','reverse')
                case 3
                    overlay_slice = squeeze(overlay_dat(:,:,slices(ix)));
                    set(handles.overlay_axes(ix),'Xdir','normal')
            end
            set(handles.overlay_axes(ix),...
                'XLim',[xmin, slice_dim(1)],...
                'YLim',[ymin, slice_dim(2)]);
            
            alpha_slice = zeros(size(overlay_slice));
            alpha_slice(overlay_slice>overlay_clim(1)) = overlay_alpha; 
            set(handles.overlay_images(ix),'CData',overlay_slice,...
                'AlphaDataMapping','none','AlphaData',alpha_slice);
        end
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY MENU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function change_background(~, ~, which_color)
    for ix = 1:10; h_background(ix).Checked = 'off'; end
    h_background(which_color).Checked = 'on';
    if which_color<10
        handles.figure.Color = background_colors(which_color,:);
        for ix = 1:length(handles.background_axes)
            set(handles.background_axes(ix),'Color',background_colors(which_color,:));
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
            set(ax_colorbar,'XLim',[0,600],'YLim',[0,100]);
        else
            figure(h_colorbar_fig);
        end
    end
    
    % Change axes color to retain contrast
    if which_color == 1 || which_color == 3 || which_color >= 5
        change_axes_color([], [], 2) 
        change_label_txt_color([], [], 2)
    elseif which_color == 2 || which_color == 4
        change_axes_color([], [], 1) 
        change_label_txt_color([], [], 1)
    end
end

function select_color(~,~,c_map)
    pt = get(gca,'CurrentPoint');
    xpos = round(pt(1,1));
    use_color = c_map(xpos,:);
    handles.figure.Color = use_color;
    for ix = 1:numslices
        set(handles.background_axes(ix),'Color',use_color);
    end
end

function change_background_thresh(~,~,~)
    % Get User Input:
    prompt = {sprintf('Voxel intensity threshold:')};
    dlg_title = 'Threshold'; num_lines = [1,35;]; 
    answer = inputdlg(prompt,dlg_title,num_lines);
    if isempty(answer); return; end
    background_thresh = str2double(answer{1});
    
    % Apply Threshold:
    for ix = 1:numslices
        background_slice = get(handles.background_images(ix),'CData');
        alpha_background = zeros(size(background_slice));
        alpha_background(background_slice>background_thresh) = 1; 
        set(handles.background_images(ix),'AlphaData',alpha_background);
    end
end

function change_show_axes(~, ~, ~)
    if show_axes
        set(h_show_axes, 'Checked', 'off'); 
        show_axes = 0;
        for ix = 1:numslices
            if isgraphics(handles.background_axes(ix),'axes')
                set(handles.background_axes(ix),'Visible','off');
            end

            % Clear TickLabels
            for ixx = 1:size(h_xticklabels, 2)
                delete(h_xticklabels(ix, ixx))
            end
            for ixx = 1:size(h_yticklabels, 2)
                delete(h_yticklabels(ix, ixx))
            end
        end
        h_xticklabels = [];
        h_yticklabels = [];
    else
        set(h_show_axes, 'Checked', 'on'); 
        show_axes = 1;
        for ix = 1:numslices
            if isgraphics(handles.background_axes(ix),'axes')
                set(handles.background_axes(ix),'Visible','on');
            end
        end
        adjust_tick_labels
    end
end

function change_x_position(~, ~, pos)
    x_tick_position = pos;
    if x_tick_position == 1 % top
        % Menu item
        set(h_x_axis_position(1), 'Checked', 'on')
        set(h_x_axis_position(2), 'Checked', 'off')
        
        % Tick labels
        x_tick_offset = 0.97;    
        
        % Axis location property
        for ix = 1:numslices
            if isgraphics(handles.background_axes(ix),'axes')
                set(handles.background_axes(ix), 'XAxisLocation','top');
            end
        end
    elseif x_tick_position == 2 % bottom
        % Menu item
        set(h_x_axis_position(1), 'Checked', 'off')
        set(h_x_axis_position(2), 'Checked', 'on')
        
        % Tick labels
        x_tick_offset = 0.04;
        
        % Axis location property
        for ix = 1:numslices
            if isgraphics(handles.background_axes(ix),'axes')
                set(handles.background_axes(ix), 'XAxisLocation','bottom');
            end
        end
    end
    
    adjust_tick_labels
end

function change_y_position(~, ~, pos)
    y_tick_position = pos;
    if y_tick_position == 1 % left
        % Menu item
        set(h_y_axis_position(1), 'Checked', 'on')
        set(h_y_axis_position(2), 'Checked', 'off')
        
        % Tick labels
        y_tick_offset = 0.02;
        
    elseif y_tick_position == 2 % right
        % Menu item
        set(h_y_axis_position(1), 'Checked', 'off')
        set(h_y_axis_position(2), 'Checked', 'on')
        
        % Tick labels
        y_tick_offset = 0.98;
        
    end
    
    adjust_tick_labels
end

function change_x_limits(~, ~, xlim)
    if nargin < 3
        
        % Obtain user input
        prompt = {sprintf('Specify X-Limits \n\n Lower:'), 'Upper:'};
        dlg_title = 'X-Limits'; num_lines = [1,20; 1,20]; 
        answer = inputdlg(prompt,dlg_title,num_lines);
        if isempty(answer); return; end
        xlim = [str2double(answer{1}), str2double(answer{2})];
        
        % Check input
        if xlim(2) < xlim(1)
           warning('Upper limit cannot be smaller than lower limit.') 
           return
        end

    elseif ~isempty(xlim)           
              
        if ~length(xlim)==2
            warning('Input ''xlim'' should be a vector of length 2.')
            return;
        end
      
    end
    
    % Convert units
    if unit_type == 2 % Physical
        xlim = xlim / slice_pixdim(1) + slice_origin(1);
    end

    % Apply limits
    for ix = 1:numslices
        set(handles.background_axes(ix),'XLim', xlim);
        if isgraphics(handles.overlay_axes(ix),'axes')
            set(handles.overlay_axes(ix),'XLim', xlim);
        end
    end
    if slice_locator_on
        set(handles.slice_locator,'YLim',xlim);
    end
    
    % Adjust ticks
    xmin = xlim(1); xmax = xlim(2);
    adjust_tick_labels
end

function change_y_limits(~, ~, ylim)
    
    if nargin < 3
        
        % Obtain user input
        prompt = {sprintf('Specify Y-Limits \n\n Lower:'),'Upper:'};
        dlg_title = 'Y-Limits'; num_lines = [1,20; 1,20]; 
        answer = inputdlg(prompt,dlg_title,num_lines);
        if isempty(answer); return; end
        ylim = [str2double(answer{1}),str2double(answer{2})];
        
        % Check input
        if ylim(2) < ylim(1)
           warning('Upper limit cannot be smaller than lower limit.') 
           return
        end
        
    elseif ~isempty(ylim)
        
        % Check input
        if ~length(ylim)==2
            warning('Input ''ylim'' should be a vector of length 2.')
            return;
        end
        
    end
    
    % Convert units
    if unit_type == 2 % Physical
        ylim = ylim / slice_pixdim(2) + slice_origin(2);
    end

    % Apply limits
    for ix = 1:numslices
        set(handles.background_axes(ix),'YLim',ylim);
        if isgraphics(handles.overlay_axes(ix),'axes')
            set(handles.overlay_axes(ix),'YLim',ylim);
        end
    end 
    if slice_locator_on
        set(handles.slice_locator,'XLim',ylim);
    end
        
    % Adjust ticks
    ymin = ylim(1); ymax = ylim(2);
    adjust_tick_labels
    
end

function change_axes_color(~,~,which_color)
    for ix = 1:3; set(h_ax_color(ix),'Checked', 'off'); end 
    set(h_ax_color(which_color), 'Checked', 'on');
    switch which_color
        case 1 % grey
            grid_color = [.5,.5,.5];
            axes_color = [.6,.6,.6];
            change_label_txt_color([], [], 1)
        case 2 % black
            grid_color = zeros(1,3);
            axes_color = zeros(1,3);
            change_label_txt_color([], [], 2)
        case 3 % white
            grid_color = ones(1,3);
            axes_color = ones(1,3);
            change_label_txt_color([], [], 3)
    end
    for ix = 1:numslices
        if isgraphics(handles.background_axes(ix),'axes')
            set(handles.background_axes(ix),'GridColor',grid_color,'MinorGridColor',grid_color,...
                'GridColorMode','manual','XColor',axes_color,'YColor',axes_color,'ZColor',axes_color);
        end
    end
    if isgraphics(handles.slice_locator,'axes')
        set(handles.slice_locator,'GridColor',grid_color,'MinorGridColor',grid_color,...
            'GridColorMode','manual','XColor',axes_color,'YColor',axes_color,'ZColor',axes_color);
    end
    % Change tick labels color
    adjust_tick_labels
end

function change_units(~, ~, unit_type_spec)
    unit_type = unit_type_spec;
    if ishandle(h_units(1)) && strcmp(get(h_units(unit_type_spec), 'Checked'),'off')
        switch unit_type
            case 1 % Voxel
                set(h_units(1), 'Checked', 'on');
                set(h_units(2), 'Checked', 'off');
            case 2 % Physical
                set(h_units(1), 'Checked', 'off');
                set(h_units(2), 'Checked', 'on');
        end
        adjust_tick_labels
        change_slice_label_position([], [], slice_label_pos)
    end
end

function change_slice_label_position(~, ~, pos)
    
    slice_label_pos = pos; 
    
    % Adjust menu checkmarks
    for ix = 1:length(h_slice_label_pos_menu)
        set(h_slice_label_pos_menu(ix),'Checked','off')
    end
    set(h_slice_label_pos_menu(slice_label_pos),'Checked','on')
   
    % Slice labels off
    if slice_label_pos == 1
        for ix = 1:numslices
            if ishandle(handles.h_slice_labels(ix))
                delete(handles.h_slice_labels(ix));
            end
        end
        return
    end
    
    % Slice labels On
    for ix = 1:numslices
        curr_ax_pos = ax_pos(ix,:);
        
        % Get position
        switch slice_label_pos
            case 2 % Top Left
                annot_pos = [curr_ax_pos(1)+.07*curr_ax_pos(3),...
                            curr_ax_pos(2)+.93*curr_ax_pos(4),...
                            .08*curr_ax_pos(3),...
                            .03];
            case 3 % Top Middle
                annot_pos = [curr_ax_pos(1)+.44*curr_ax_pos(3),...
                            curr_ax_pos(2)+.93*curr_ax_pos(4),...
                            .14*curr_ax_pos(3),...
                            .03];
            case 4 % Top Right
                annot_pos = [curr_ax_pos(1)+.86*curr_ax_pos(3),...
                            curr_ax_pos(2)+.93*curr_ax_pos(4),...
                            .08*curr_ax_pos(3),...
                            .03]; 
            case 5 % Bottom Left
                annot_pos = [curr_ax_pos(1)+.07*curr_ax_pos(3),...
                            curr_ax_pos(2),...
                            .08*curr_ax_pos(3),...
                            .03];
            case 6 % Bottom Middle
                annot_pos = [curr_ax_pos(1)+.44*curr_ax_pos(3),...
                            curr_ax_pos(2),...
                            .14*curr_ax_pos(3),...
                            .04];
            case 7 % Bottom Right
                annot_pos = [curr_ax_pos(1)+.86*curr_ax_pos(3),...
                            curr_ax_pos(2),...
                            .08*curr_ax_pos(3),...
                            .03];
        end
        
        % Get label
        if unit_type == 2
            slice_distances(ix) = round((slices(ix)-slice_origin(3))*slice_pixdim(3),3,'significant');
            curr_str = sprintf('%2.1f',round(slice_distances(ix),1));
        else
            curr_str = num2str(slices(ix));
        end
        
        if ~ishandle(handles.h_slice_labels(ix))% nargin>3 && initial_call
            handles.h_slice_labels(ix) = annotation('textbox','Position',annot_pos,...
                'String',curr_str,'FontName','Helvetica','FontSize',14,...
                'EdgeColor','none','HorizontalAlignment','center',...
                'VerticalAlignment','middle', 'color', slice_label_text_color,...
                'FitBoxToText','on','Margin',0);
        else
            set(handles.h_slice_labels(ix),'Position',annot_pos,...
                'String',curr_str,'FontSize',14,'HorizontalAlignment','center',...
                'FitBoxToText','on','Margin',0, 'color', slice_label_text_color)
        end
        
    end
end

function change_slice_labels_background(~, ~, which_color)
   for ix = 1:9; h_slice_label_background_color(ix).Checked = 'off'; end
    h_slice_label_background_color(which_color).Checked = 'on';
    if which_color==1
        for ix = 1:numslices
            if ishandle(handles.h_slice_labels(ix))
                set(handles.h_slice_labels(ix),'BackgroundColor','none')
            end
        end
    elseif which_color < 10
        cmap = annotation_background_colors(which_color,:);
        for ix = 1:numslices
            if ishandle(handles.h_slice_labels(ix))
                set(handles.h_slice_labels(ix),'BackgroundColor',cmap)
            end
        end
    else
        if ~ishandle(h_overlay_colorbar_fig)
            h_overlay_colorbar_fig = figure('menu','none','Name','Choose Color','NumberTitle','off');
            h_overlay_colorbar_fig.Position(4) = 70; 
            ax_colorbar = gca; set(ax_colorbar,'Position',[0,0,1,1],'Visible','off','nextplot','add');
            c_map = jet(600); c_image = zeros(100,600,3);
            c_image(:,:,1) = repmat(c_map(:,1)',100,1);
            c_image(:,:,2) = repmat(c_map(:,2)',100,1);
            c_image(:,:,3) = repmat(c_map(:,3)',100,1);
            image(c_image,'Parent',ax_colorbar,'ButtonDownFcn',...
                {@manual_select_labels_background_color,c_map});
            set(ax_colorbar,'XLim',[0,600],'YLim',[0,100])
        else
            figure(h_overlay_colorbar_fig)
        end
    end
end

function manual_select_labels_background_color(~,~,c_map)
    pt = get(gca,'CurrentPoint');
    xpos = round(pt(1,1));
    use_color = c_map(xpos,:);
    for ix = 1:numslices
        if ishandle(handles.h_slice_labels(ix))
            set(handles.h_slice_labels(ix),'BackgroundColor',use_color)
        end
    end
end

function change_label_txt_color(~,~,which_color)
    for ix = 1:3; set(h_slice_label_color(ix), 'Checked','off'); end
    set(h_slice_label_color(which_color), 'Checked','on')
    switch which_color
        case 1 % grey 
            slice_label_text_color = [.5,.5,.5];
            for ix = 1:numslices
                if ishandle(handles.h_slice_labels(ix))
                    set(handles.h_slice_labels(ix),'Color',slice_label_text_color); 
                end
            end
        case 2 % black
            slice_label_text_color = zeros(1,3);
            for ix = 1:numslices
                if ishandle(handles.h_slice_labels(ix))
                    set(handles.h_slice_labels(ix),'Color',slice_label_text_color); 
                end
            end
        case 3 % white
            slice_label_text_color = ones(1,3);
            for ix = 1:numslices
                if ishandle(handles.h_slice_labels(ix))
                    set(handles.h_slice_labels(ix),'Color',slice_label_text_color); 
                end
            end
    end
end

function adjust_tick_labels

    if ~show_axes; return; end
    
    % Check measurements
    % (slice_dim - slice_origin) .* slice_pixdim
    % ([1,1,1] - slice_origin) .* slice_pixdim
    
    axes_properties = get(handles.background_axes(1));

    xticks = round(linspace(xmin, xmax, n_ticks_x + 2));
    yticks = round(linspace(ymin, ymax, n_ticks_y + 2)); 
    
    xticks = xticks(2:end-1);
    yticks = yticks(2:end-1);
    
    xtick_norm = (xticks - axes_properties.XLim(1)) / (axes_properties.XLim(2) - axes_properties.XLim(1));
    ytick_norm = (yticks - axes_properties.YLim(1)) / (axes_properties.YLim(2) - axes_properties.YLim(1));
    
    if unit_type == 1 % Voxel
        xticklabels = sprintfc('%1g', round(xticks));
        yticklabels = sprintfc('%1g', round(yticks));
    elseif unit_type == 2 % Physical
        xticklabels = sprintfc('%.1f', round((xticks - slice_origin(1)) .* slice_pixdim(1), 3, 'significant'));
        yticklabels = sprintfc('%.1f', round((yticks - slice_origin(2)) .* slice_pixdim(2), 3, 'significant'));
    end
    
    if dimension == 2 % Sagittal needs x-reversed
       xticklabels = fliplr(xticklabels); 
    end
    
    for ix = 1:length(handles.background_axes)
        
        set(handles.background_axes(ix), ...
            'XTick', xticks, ...
            'YTick', yticks, ...
            'XTickLabel', {}, ...
            'YTickLabel', {});
        
        if x_tick_position == 1
            set(handles.background_axes(ix), ...
                'XAxisLocation', 'top')
        elseif x_tick_position == 2
            set(handles.background_axes(ix), ...
                'XAxisLocation', 'bottom')
        end
        
        if y_tick_position == 1
            set(handles.background_axes(ix), ...
                'YAxisLocation', 'left')
        elseif y_tick_position == 2
            set(handles.background_axes(ix), ...
                'YAxisLocation', 'right')
        end
        
        if size(h_xticklabels, 1) >= ix
            if isgraphics(h_xticklabels(ix, 1), 'text')
                for ixx = 1:size(h_xticklabels,2)
                    set(h_xticklabels(ix, ixx), ...
                        'Position', [xtick_norm(ixx), x_tick_offset], ...
                        'String', xticklabels(ixx), ... 
                        'units', 'normalized', ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', axes_properties.FontSize, ...
                        'FontName', axes_properties.FontName, ...
                        'FontWeight', axes_properties.FontWeight, ...
                        'FontAngle', axes_properties.FontAngle, ...
                        'color', axes_properties.XColor, ...
                        'parent', handles.background_axes(ix));
                end
            end
        else 
            h_xticklabels(ix, :) = text(xtick_norm, ...
                repmat(x_tick_offset, 1, n_ticks_x), xticklabels,... 
                'units', 'normalized', ...
                'HorizontalAlignment', 'center', ...
                'FontSize', axes_properties.FontSize, ...
                'FontName', axes_properties.FontName, ...
                'FontWeight', axes_properties.FontWeight, ...
                'FontAngle', axes_properties.FontAngle, ...
                'color', axes_properties.XColor, ...
                'parent', handles.background_axes(ix));
        end

        if size(h_yticklabels, 1) >= ix
            if isgraphics(h_yticklabels(ix, 1), 'text')
                for ixx = 1:size(h_yticklabels, 2)
                    set(h_yticklabels(ix, ixx), ...
                        'Position', [y_tick_offset, ytick_norm(ixx)], ...
                        'String', yticklabels(ixx),...
                        'units', 'normalized', ...
                        'VerticalAlignment', 'middle', ...
                        'FontSize', axes_properties.FontSize, ...
                        'FontName', axes_properties.FontName, ...
                        'FontWeight', axes_properties.FontWeight, ...
                        'FontAngle', axes_properties.FontAngle, ...
                        'color', axes_properties.XColor, ...
                        'parent', handles.background_axes(ix));
                end
            end
        else 
            h_yticklabels(ix, :) = text(repmat(y_tick_offset, 1, n_ticks_y), ...
                ytick_norm, yticklabels,...
                'units', 'normalized', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', axes_properties.FontSize, ...
                'FontName', axes_properties.FontName, ...
                'FontWeight', axes_properties.FontWeight, ...
                'FontAngle', axes_properties.FontAngle, ...
                'color', axes_properties.XColor, ...
                'parent', handles.background_axes(ix));
        end
        
        % Adjust text alignment based on axes position:
        if size(h_xticklabels, 1) >= ix
            if isgraphics(h_xticklabels(ix, 1), 'text')
                for ixx = 1:size(h_xticklabels,2)
                    if x_tick_position == 1
                        set(h_yticklabels(ix, ixx), ...
                            'VerticalAlignment', 'top')
                    elseif x_tick_position == 2
                        set(h_yticklabels(ix, ixx), ...
                            'VerticalAlignment', 'bottom')
                    end
                end
            end
        end
        
        if size(h_yticklabels, 1) >= ix
            if isgraphics(h_yticklabels(ix, 1), 'text')
                for ixx = 1:size(h_yticklabels, 2)
                    if y_tick_position == 1
                        set(h_yticklabels(ix, ixx), ...
                            'VerticalAlignment', 'middle', ...
                            'HorizontalAlignment', 'left')
                    elseif y_tick_position == 2
                        set(h_yticklabels(ix, ixx), ...
                            'VerticalAlignment', 'middle', ...
                            'HorizontalAlignment', 'right')
                    end
                end
            end
        end

    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OVERLAY MENU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function change_overlay_clim(~,~,from_other)
    if nargin<3 || isempty(from_other)
        prompt = {sprintf('Specify color limits: \n\nLower Bound:'),'Upper Bound:'};
        dlg_title = 'CLim'; num_lines = [1,35;1,35]; 
        answer = inputdlg(prompt,dlg_title,num_lines);
        if isempty(answer); return; end
        overlay_clim(1) = str2double(answer{1}); overlay_clim(2) = str2double(answer{2});
        for ix = 1:numslices
            overlay_slice = overlay_slice_data(:,:,ix);
            background_slice= background_slice_data(:,:,ix);
            alpha_slice = zeros(size(overlay_slice));
            alpha_slice(overlay_slice>overlay_clim(1)) = overlay_alpha;
%             alpha_slice(overlay_slice>overlay_clim(2)) = 0;
            alpha_slice(background_slice<=background_thresh) = 0;
            overlay_vec = overlay_slice(:);
%             ind_nan = ((overlay_vec<=overlay_clim(1))+(overlay_vec>overlay_clim(2)))>0;
            ind_nan = overlay_vec<=overlay_clim(1);
            overlay_vec = min(m,round((m-1).*(overlay_vec-overlay_clim(1))./(overlay_clim(2)-overlay_clim(1)))+1);
            overlay_vec(overlay_vec<=0) = 1; % assure no negative or 0 indices
            overlay_vec(ind_nan) = nan;
            overlay_slice(:) = overlay_vec;
            set(handles.overlay_images(ix),'CData',overlay_slice,'AlphaData',alpha_slice,'CDataMapping','direct')
        end
    end
end
    
function overlay_opacity_callback(~, ~, which_alpha)
    switch which_alpha
        case 1; overlay_alpha = 1;
        case 2; overlay_alpha = .9;
        case 3; overlay_alpha = .8;
        case 4; overlay_alpha = .7;    
        case 5; overlay_alpha = .6;
        case 6; overlay_alpha = .5;
        case 7; overlay_alpha = .4;
        case 8; overlay_alpha = .3;
        case 9; overlay_alpha = .2;
        case 10; overlay_alpha = .1;
        case 11; overlay_alpha = 0;
    end
    for ix = 1:numslices
        alpha_overlay = get(handles.overlay_images(ix),'AlphaData');
        alpha_overlay(alpha_overlay~=0) = overlay_alpha;
        set(handles.overlay_images(ix),'AlphaData',alpha_overlay);
    end
end

function change_overlay_single_colors(~,~,which_color)
    for ix = 1:8; h_overlay_single_colors(ix).Checked = 'off'; end
    h_overlay_single_colors(which_color).Checked = 'on';
    if which_color < 8
        cmap = roi_single_colors(which_color,:);
        for ix = 1:numslices
            colormap(handles.overlay_axes(ix),cmap)
        end
    else
        if ~ishandle(h_overlay_colorbar_fig)
            h_overlay_colorbar_fig = figure('menu','none','Name','Choose Color','NumberTitle','off');
            h_overlay_colorbar_fig.Position(4) = 70; 
            ax_colorbar = gca; set(ax_colorbar,'Position',[0,0,1,1],'Visible','off','nextplot','add');
            c_map = jet(600); c_image = zeros(100,600,3);
            c_image(:,:,1) = repmat(c_map(:,1)',100,1);
            c_image(:,:,2) = repmat(c_map(:,2)',100,1);
            c_image(:,:,3) = repmat(c_map(:,3)',100,1);
            image(c_image,'Parent',ax_colorbar,'ButtonDownFcn',{@manual_select_overlay_color,c_map});
            set(ax_colorbar,'XLim',[0,600],'YLim',[0,100])
        else
            figure(h_overlay_colorbar_fig)
        end
    end
end

function manual_select_overlay_color(~,~,c_map)
    pt = get(gca,'CurrentPoint');
    xpos = round(pt(1,1));
    use_color = c_map(xpos,:);
    for ix = 1:numslices
        colormap(handles.overlay_axes(ix),use_color)
    end
end

function overlay_colormap_callback(hObject,~,which_color)
    for ix = 1:14; set(h_overlay_color_spec(ix), 'Checked', 'off'); end
    set(h_overlay_color_spec(which_color), 'Checked', 'on');
    if ishandle(hObject)
        cmap_name = get(hObject,'Label');
    elseif ischar(hObject)
        cmap_name = hObject;    
    end
    switch cmap_name
        case 'Blue-White-Red'
            cmap = bluewhitered(m,overlay_clim(1),overlay_clim(2));
        case 'Red-Blue'
            cmap = redblue(m);
        otherwise
        cmap = eval([cmap_name,'(',num2str(m),')']);
    end
    for ix = 1:numslices
        colormap(handles.overlay_axes(ix),cmap)
    end
    % Change back to scaled colormapping if roi_colors was specified
    if ~isempty(roi_colors) && first_cmap_change
        for ix = 1:numslices
%             overlay_slice = get(handles.overlay_images(ix),'CData');
            overlay_slice = overlay_slice_data(:,:,ix);
            m = 200; 
            overlay_vec = overlay_slice(:);
            ind_nan = ((overlay_vec<=overlay_clim(1))+(overlay_vec>overlay_clim(2)))>0;
            overlay_vec = min(m,round((m-1).*(overlay_vec-overlay_clim(1))./(overlay_clim(2)-overlay_clim(1)))+1);
            overlay_vec(overlay_vec<=0) = 1; % assure no negative or 0 indices
            overlay_vec(ind_nan) = nan; 
            overlay_slice(:) = overlay_vec;
            set(handles.overlay_images(ix),'CData',overlay_slice)
        end
        first_cmap_change = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLICE LOCATOR CALLBACKS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mesh_face_opacity_callback(~, ~, which_alpha)
    for ix = 1:11; set(mesh_face_alpha_submenus(ix),'Checked','off'); end
    set(mesh_face_alpha_submenus(which_alpha),'Checked','on') 
    switch which_alpha
        case 1; mesh_alpha = 1;
        case 2; mesh_alpha = .9;
        case 3; mesh_alpha = .8;
        case 4; mesh_alpha = .7;    
        case 5; mesh_alpha = .6;
        case 6; mesh_alpha = .5;
        case 7; mesh_alpha = .4;
        case 8; mesh_alpha = .3;
        case 9; mesh_alpha = .2;
        case 10; mesh_alpha = .1;
        case 11; mesh_alpha = 0;
    end
    for ix = 1:numslices
        set(handles.mesh(ix),'AlphaData',mesh_alpha,'FaceAlpha',mesh_alpha);
    end
end

function mesh_edge_opacity_callback(~, ~, which_alpha)
    for ix = 1:11; set(mesh_edge_alpha_submenus(ix),'Checked','off'); end
    set(mesh_edge_alpha_submenus(which_alpha),'Checked','on') 
    switch which_alpha
        case 1; mesh_alpha = 1;
        case 2; mesh_alpha = .9;
        case 3; mesh_alpha = .8;
        case 4; mesh_alpha = .7;    
        case 5; mesh_alpha = .6;
        case 6; mesh_alpha = .5;
        case 7; mesh_alpha = .4;
        case 8; mesh_alpha = .3;
        case 9; mesh_alpha = .2;
        case 10; mesh_alpha = .1;
        case 11; mesh_alpha = 0;
    end
    for ix = 1:numslices
        set(handles.mesh(ix),'EdgeAlpha',mesh_alpha);
    end
end

%% ADDITIONAL UTILITIES:

function resizeFigure
    aspect_ratio = (slice_dim(1) * slice_pixdim(1) / ax_pos(3)) / (slice_dim(2) * slice_pixdim(2) / ax_pos(4));
    if aspect_ratio > 1
        set(handles.figure, 'position', [0, .5*(1-(1/aspect_ratio)), 1, 1/aspect_ratio])
    else
        set(handles.figure, 'position', [.5*(1-aspect_ratio), 0, aspect_ratio, 1])
    end
end

function status = load_new_background
    status = 0;
    [background, background_path] = uigetfile({'*.nii';'*.nii.gz';'*.img'},...
        'Select a 3D image:','MultiSelect','off');
    if background_path==0 % Cancel
        disp('Load new image cancelled.'); 
        return
    end
    if ischar(background)
        try 
            back_img = load_nii(fullfile(background_path,background)); 
        catch
            try
                back_img = load_untouch_nii(fullfile(background_path,background));
                warning('Non-orthogonal shearing detected in affine matrix of background image. Loaded successfully without applying affine.')
            catch
                error('Error:  failed to load background image')
            end
        end
        imdat = back_img.img; 
        pixdim = back_img.hdr.dime.pixdim([2,3,4]);
        origin = back_img.hdr.originator;
        switch bitand(back_img.hdr.dime.xyzt_units, 7) % see xform_nii.m, extra_nii_hdr.m
            case 1, physical_units = '(m)';
            case 2, physical_units = '(mm)';
            case 3, physical_units = '(microns)';
            otherwise, physical_units = '';
        end
    end
    status = 1;
end

function h = neuroimage_overobj(Type)
    %   Adapted from internal Matlab function overobj.m
    %   H = OVEROBJ(TYPE) check searches visible objects of Type TYPE in 
    %   the PointerWindow looking for one that is under the pointer.  It
    %   returns the handle to the first object it finds under the pointer
    %   or else the empty matrix.
    %   Copyright 1984-2013 The MathWorks, Inc.

    fig = matlab.ui.internal.getPointerWindow();
    % Look for quick exit
    if fig==0
       h = [];
       return
    end

    % Assume root units are pixels
    p = get(0,'PointerLocation');
    % Get figure position in pixels
    figUnit = get(fig,'Units');
    set(fig,'Units','pixels');
    figPos = get(fig,'Position');
    set(fig,'Units',figUnit)

    x = (p(1)-figPos(1))/figPos(3);
    y = (p(2)-figPos(2))/figPos(4);
    c = findobj(get(fig,'Children'),'flat','Type',Type); % previously only searched visible objects
    for h = c'
       hUnit = get(h,'Units');
       set(h,'Units','norm')
       r = get(h,'Position');
       set(h,'Units',hUnit)
       if ( (x>r(1)) && (x<r(1)+r(3)) && (y>r(2)) && (y<r(2)+r(4)) )
          return
       end
    end
    h = [];
end

end % end nifti_studio_mosaic.m