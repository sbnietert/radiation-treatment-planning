function voxels = setSurfaceEnclosedVolume(voxels, coordinates, ...
                                           surface, value)
%SETSURFACEENCLOSEDVOLUME Set region enclosed by surface to value.
%   VOXELS = SETSURFACEENCLOSEDVOLUME(VOXELS, COORDINATES, SURFACE, VALUE)
%   returns voxel grid described by COORDINATES with volume enclosed by
%   SURFACE set to VALUE.
%
%   See also COORDINATEDATA.

if coordinates.Dim == 2
    iSurface = surface(:, 1);
    jSurface = surface(:, 2);
    iVoxels = coordinates.I(:);
    jVoxels = coordinates.J(:);
    inSurface = inpolygon(iVoxels, jVoxels, iSurface, jSurface);
    voxels(inSurface) = value;
else
    error('setSurfaceEnclosedVolume only supports 2D grids')
end
end