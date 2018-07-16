function voxels = setAABB(voxels, coordinates, corner1, corner2, value)
%SETAABB Set all voxels within axis-aligned bounding box to value.
%   VOXELS = SETAABB(VOXELS, COORDINATES, CORNER1, CORNER2, VALUE) returns
%   voxel grid described by COORDINATES with all voxels between lower
%   corner CORNER1 and upper corner CORNER2 to VALUE.
%
%   See also COORDINATEDATA.

if or(length(corner1) ~= coordinates.Dim, ...
      length(corner2) ~= coordinates.Dim)
    error('corner and coordinate dimensions do not match');
end

%- Create logical filter for box
filter1 = and(coordinates.I > corner1(1), ...
              coordinates.I < corner2(1));
filter2 = and(coordinates.J > corner1(2), ...
              coordinates.J < corner2(2));
if coordinates.Dim == 3
    filter3 = and(coordinates.K > corner1(3), ...
              coordinates.K < corner2(3));
else
    filter3 = ones(coordinates.GridSize);
end
filter = and(filter1, and(filter2, filter3));

%- Set all voxels within box to value
voxels(filter) = value;
end