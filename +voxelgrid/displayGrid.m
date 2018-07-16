function displayGrid(voxels)
%DISPLAYGRID Display voxel grid.
%   DISPLAYGRID(VOXELS) displays (slice of) VOXELS using imagesc, with
%   custom data cursor and slice slider, if grid is 3D.

%- Select data and set up axes
gridSize = size(voxels);
dim = length(gridSize);
f = figure;
if dim == 2
    data = voxels;
    a = axes(f);
else
    k = round(gridSize(3) / 2);
    data = squeeze(voxels(:, :, k));
    a = axes(f, 'Position', [0.05 0.2 1 0.75]);
end

%- Display data
function drawSlice
    imagesc(data);
    set(a, 'xtick', []);
    set(a, 'ytick', []);
    pbaspect(a, [1, 1, 1]);
    caxis([0 max(voxels(:))]);
    colormap default;
    colorbar;
end
drawSlice()

%- Set up data cursor
function txt = dataCursorUpdate(~, event)
    data = get(event, 'Position');
    i = data(2);
    j = data(1);
    if dim == 2
        txt = sprintf('(%d, %d)', i, j);
    else
        txt = sprintf('(%d, %d, %d)', i, j, k);
    end
end
dcm = datacursormode(f);
set(dcm, 'Enable', 'on');
set(dcm, 'UpdateFcn', @dataCursorUpdate);

%- No slider for 2D grid
if dim == 2
    return
end

%- Set up slider for 3D grid
function sliderUpdate(handle, ~)
    k = round(handle.Value);
    data = squeeze(voxels(:, :, k));
    drawSlice()
end
slider = uicontrol('Parent', f, 'Style', 'slider', ...
                    'Units', 'Normalized', ...
                    'Position', [0.25 0.05 .5 0.1], ...
                    'value', k, ...
                    'min', 1, 'max', gridSize(3));
slider.Callback = @sliderUpdate;
end