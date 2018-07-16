function isocenters = genIsocenters(targetVolume, coordinates, ...
                                    numIsocenters, collimatorSizes, ...
                                    alpha, f)
%GENISOCENTERS Generate isocenters for target volume.
%   ISOCENTERS = GENISOCENTERS(TARGETVOLUME, COORDINATES, NUMISOCENTERS,
%   COLLIMATORSIZES, ALPHA, F) returns up to NUMISOCENTERS isocenter
%   coordinates generated for volume defined by binary voxel mask
%   TARGETVOLUME with coordinate data COORDINATES using collimator sizes
%   COLLIMATORSIZES. Closeness acceptance factor ALPHA and relaxation
%   factor F are hyperparameters which control the scoring/selection
%   process described in https://www.ncbi.nlm.nih.gov/pubmed/22755698.
%
%   See also GETGRASSFIREDISTANCES, GETCANDIDATESCORES.

isocenters = [];
isocenterSizes = [];
currentTarget = targetVolume;

for i = 1:numIsocenters
    %- Determine next set of isocenter candidates
    G = planning.getGrassfireDistances(currentTarget);
    %{
    voxelgrid.displayGrid(G);
    pause; close all;
    %}
    delta = max(G(:)); % maximum grassfire depth
    candidateIndices = find(G == delta); % indices of deepest voxels
    candidatesCell = cell(coordinates.Dim, 1);
    [candidatesCell{:}] = ind2sub(size(G), candidateIndices);
    candidates = cell2mat(candidatesCell'); % candidate voxel coordinates
    numCandidates = length(candidateIndices); % number of candidates
    
    %- Compute max collimator size for each candidate
    candidateSizes = zeros(numCandidates, 1);
    for j = 1:numCandidates
        sizeFilter = (1 - f) * collimatorSizes <= ...
                     2 * G(candidateIndices(j));
        selectedSize = max(collimatorSizes(sizeFilter));
        if isempty(selectedSize)
            % => all collimator sizes are too large for candidate
            candidateSizes(j) = -1;
        else
            candidateSizes(j) = selectedSize;
        end
    end
    
    %- Only pursue candidates with valid collimator size
    candidateFilter = candidateSizes >= 0;
    candidates = candidates(candidateFilter, :);
    candidateSizes = candidateSizes(candidateFilter);
    if isempty(candidates)
        break;
    end
    
    %- Compute score for each candidate
    scores = planning.getCandidateScores(candidates, candidateSizes, ...
                                         isocenters, isocenterSizes, ...
                                         currentTarget, delta, alpha, f);
    
    %- Select best candidate to be isocenter
    [~, argmax] = max(scores);
    selectedIsocenter = candidates(argmax, :);
    selectedSize = candidateSizes(argmax);
    isocenters(end+1, :) = selectedIsocenter;
    isocenterSizes(end+1) = selectedSize;
    
    %- Remove covered voxels from current target volume
    currentTarget = voxelgrid.setSphere(currentTarget, coordinates, ...
                                        selectedIsocenter, ...
                                        selectedSize, 0);
end
end

