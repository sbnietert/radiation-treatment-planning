function scores = getCandidateScores(candidates, candidateSizes, ...
                                     isocenters, isocenterSizes, ...
                                     targetVolume, delta, alpha, f)
%GETCANDIDATESCORES Compute scores for isocenter candidates.
%   SCORES = GETCANDIDATESCORES(CANDIDATES, CANDIDATESIZES, ISOCENTERS,
%   ISOCENTERSIZES, TARGETVOLUME, DELTA, ALPHA, F) returns the score for
%   each candidate i with coordinates CANDIDATES(i) and size
%   CANDIDATESIZES(i), given previously selected isocenter coordinates
%   ISOCENTERS and sizes ISOCENTERSIZES, within volume defined by binary
%   voxel mask TARGETVOLUME with grassfire depth DELTA. Closeness
%   acceptance factor ALPHA and relaxation factor F are hyperparameters
%   which control the scoring process described in
%   https://www.ncbi.nlm.nih.gov/pubmed/22755698.
%
%   See also GENISOCENTERS.

numCandidates = size(candidates, 1);
numIsocenters = size(isocenters, 1);
gridSize = size(targetVolume);
dim = length(gridSize);

%- Compute boundary scores
depthScores = zeros(numCandidates, 1);
for i = 1:numCandidates
    candidate = num2cell(candidates(i, :));
    
    %- Compute distance to boundary in "smaller" directions
    smallDists = zeros(dim, 1);
    for d = 1:dim
        smaller = candidate;
        for s = candidate{d}:-1:1
            smaller{d} = s;
            if ~targetVolume(sub2ind(gridSize, smaller{:}))
                break
            end
        end
        smallDists(d) = candidate{d} - s;
    end
    
    %- Compute distance to boundary in "larger" directions
    largeDists = zeros(dim, 1);
    for d = 1:dim
        larger = candidate;
        for l = candidate{d}:-1:1
            larger{d} = l;
            if ~targetVolume(sub2ind(gridSize, larger{:}))
                break
            end
        end
        largeDists(d) = l - candidate{d};
    end
    
    depths = max(smallDists, largeDists);
    depthScores(i) = mean(depths) / delta;
end
boundaryScores = depthScores / max(depthScores);

%- Compute iso scores if previous isocenters have been selected
if numIsocenters == 0
    scores = boundaryScores;
else
    %- Compute pairwise distances between candidates and isocenters
    pairwiseDistances = zeros(numCandidates, numIsocenters);
    for i = 1:numCandidates
        candidate = candidates(i, :);
        offsets = isocenters - candidate;
        pairwiseDistances(i, :) = sqrt(sum(offsets.^2, 2));
    end

    %- Compute iso scores
    pairwiseIsoScores = zeros(numCandidates, numIsocenters);
    for i = 1:numCandidates
        for k = 1:numIsocenters
            dist = pairwiseDistances(i, k);
            maxDist = max(pairwiseDistances(i, :));
            optDist = (1 + f) * (candidateSizes(i) + isocenterSizes(k));
            if dist < alpha * optDist
                pairwiseIsoScores(i, k) = 0;
            else
                pairwiseIsoScores(i, k) = dist / maxDist;
            end
        end
    end
    isoScores = min(pairwiseIsoScores, [], 2);
    scores = boundaryScores .* isoScores;
end
end