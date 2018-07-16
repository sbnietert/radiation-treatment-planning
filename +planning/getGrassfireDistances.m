function G = getGrassfireDistances(targetVolume)
%GETGRASSFIREDISTANCES Compute grassfire transform of target volume.
%   G = GETGRASSFIREDISTANCES(TARGETVOLUME) returns the grassfire transform
%   (a.k.a. L1 distance transform) of the volume defined by binary voxel
%   mask TARGETVOLUME, with G(v) = L1 distance of voxel v to boundary of
%   volume for v in volume, 0 otherwise.
%
%   See the following page for more details:
%   https://en.wikipedia.org/wiki/Grassfire_transform.
%
%   See also GENISOCENTERS.

numVoxels = numel(targetVolume);
gridSize = size(targetVolume);
dim = length(gridSize);
G = zeros(gridSize);

%- Perform forward sweep through voxels, setting G for target voxels v to
%  G(v) = 1 + min_u G(u), minimizing over previously visited neighbors u
for v = 1:numVoxels
    if targetVolume(v)
        coord = cell(dim, 1);
        [coord{:}] = ind2sub(gridSize, v);
        
        n = [];
        for d = 1:dim
            if coord{d} > 1
                prev = coord;
                prev{d} = prev{d} - 1;
                n(end+1) = G(sub2ind(gridSize, prev{:}));
            end
        end

        G(v) = 1;
        if ~isempty(n)
            G(v) = G(v) + min(n);
        end
    end
end

%- Perform backwards sweep through voxels, setting G for target voxels v to
%  G(v) = min(G(v), 1 + min_u G(u)), minimizing over previously visited
%  neighbors u
for v = numVoxels:-1:1
    if targetVolume(v)
        coord = cell(dim, 1);
        [coord{:}] = ind2sub(gridSize, v);
        
        n = [];
        for d = 1:dim
            if coord{d} < gridSize(d)
                next = coord;
                next{d} = next{d} + 1;
                n(end+1) = G(sub2ind(gridSize, next{:}));
            end
        end

        if ~isempty(n)
            G(v) = min(G(v), 1 + min(n));
        end
    end
end
end

