function doseRate = getDoseRate(coordinates, isocenter, collimatorSize)
%GETDOSERATE Compute dose rate from shot for all voxels in grid.
%   DOSERATE = GETDOSERATE(COORDINATES, ISOCENTER, COLLIMATORSIZE) returns
%   the dose rate from radiation shot of size COLLIMATORSIZE administered
%   at ISOCENTER for all voxels in grid described by COORDINATES. The dose
%   rate is assumed to be a Gaussian with mean ISOCENTER and standard 
%   deviation COLLIMATORSIZE.

    %- Compute squared distance of each voxel to isocenter
    dI = coordinates.I - isocenter(1);
    dJ = coordinates.J - isocenter(2);
    D2 = dI.^2 + dJ.^2;
    if coordinates.Dim == 3
        dK = coordinates.K - isocenter(3);
        D2 = D2 + dK.^2;
    end
    
    %- Evaluate Gaussian
    doseRate = 1 / sqrt(2 * pi * collimatorSize^2) * ...
               exp(-D2 / (2 * collimatorSize^2));
end