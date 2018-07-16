function voxels = setSphere(voxels, coordinates, center, radius, value)
%SETSPHERE Set sphere of voxels to value.
%   VOXELS = SETSPHERE(VOXELS, COORDINATES, CENTER, RADIUS,
%   VALUE) returns voxel grid described by COORDINATES with
%   sphere of radius RADIUS and center CENTER set to VALUE.
%
%   See also COORDINATEDATA.

if length(center) ~= coordinates.Dim
    error('center and coordinate dimensions do not match');
end

%- Compute squared distance of each voxel to center
dI = coordinates.I - center(1);
dJ = coordinates.J - center(2);
D2 = dI.^2 + dJ.^2;
if coordinates.Dim > 2
    dK = coordinates.K - center(3);
    D2 = D2 + dK.^2;
end

%- Set all voxels within radius of center to value
voxels(D2 <= radius*radius) = value;
end