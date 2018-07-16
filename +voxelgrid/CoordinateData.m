classdef CoordinateData
    %CoordinateData Stores coordinates for voxel grid
    %   This data is useful for performing computations on all voxels
    %   coordinates at the same time.
    
    properties (GetAccess = public, SetAccess = private)
        I % i-coordinates, i.e. I(i,j)/I(i,j,k) = i
        J % j-coordinates
        K % k-coordinates, only for 3D grids
        GridSize % size of voxel grid
        Dim % length(GridSize) = 2 or 3
    end
    
    methods
        function obj = CoordinateData(gridSize)
            % Create coordinate data for voxel grid of specified size
            obj.GridSize = gridSize;
            obj.Dim = length(gridSize);
            if obj.Dim == 2
                [obj.J, obj.I] = meshgrid(1:gridSize(2), ...
                                          1:gridSize(1));
            elseif obj.Dim == 3
                [obj.J, obj.I, obj.K] = meshgrid(1:gridSize(2), ...
                                                 1:gridSize(1), ...
                                                 1:gridSize(3));
            else
                error('only 2D and 3D voxel grids supported');
            end
        end
    end
end

