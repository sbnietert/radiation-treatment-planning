function plan = genTreatmentPlan(structureAssignments, ...
                                                 minDosePerStructure, ...
                                                 maxDosePerStructure, ...
                                                 underdoseWeights, ...
                                                 overdoseWeights, ...
                                                 doseRates)
%GENTREATMENTPLAN Compute optimal duration of each shot.
%   PLAN = GENTREATMENTPLAN(STRUCTUREASSIGNMENTS, MINDOSEPERSTRUCTURE,
%   MAXDOSEPERSTRUCTURE, UNDERDOSEWEIGHTS, OVERDOSEWEIGHTS, DOSERATES) 
%   returns the treatment plan with optimal duration for each shot
%   administered during radiation treatment for voxel grid with structures
%   given by STRUCTUREASSIGNMENTS, structure dose limits given by
%   MINDOSEPERSTRUCTURE and MAXDOSEPERSTRUCTURE, extreme dose weightings
%   given by UNDERDOSEWEIGHTS and OVERDOSEWEIGHTS, and dose rates for each
%   shot given by DOSERATES. Durations computed using quadratic program
%   outlined in https://tspace.library.utoronto.ca/handle/1807/34010.
    
    %- Get min and max dose for each voxel, compute structure voxel counts
    gridSize = size(structureAssignments);
    minDoses = zeros(gridSize);
    maxDoses = zeros(gridSize);
    numStructures = length(minDosePerStructure);
    structureVoxelCounts = zeros(numStructures, 1);
    for s = 1:numStructures
        structureFilter = structureAssignments == s;
        minDoses(structureFilter) = minDosePerStructure(s);
        maxDoses(structureFilter) = maxDosePerStructure(s);
        structureVoxelCounts(s) = sum(structureFilter(:));
    end
    
    %- Find active voxels and get counts of various objects
    activeVoxels = find(structureAssignments > 0);
    numActiveVoxels = length(activeVoxels);
    numIsocenters = size(doseRates, 1);
    numCollimatorSizes = size(doseRates, 2);
    numShots = numIsocenters * numCollimatorSizes;
    
    %- See paper for description of objective function and constraints
    %- 2 variables per active voxel, 1 variable per shot
    H = speye(2*numActiveVoxels + numShots);
    H((2*numActiveVoxels + 1):end,(2*numActiveVoxels + 1):end) = 0;
    f = zeros(numActiveVoxels*2 + numShots, 1);
    
    %- 1 underdose constraint and 1 overdose constraint per voxel
    b = zeros(2*numActiveVoxels, 1);
    
    %- Precompute entries of sparse matrix A
    rows = zeros(2*numActiveVoxels*(1 + numShots), 1);
    cols = zeros(size(rows));
    vals = zeros(size(rows));
    entryIndex = 1;
    for v = 1:numActiveVoxels
        voxelIndex = activeVoxels(v);
        structure = structureAssignments(voxelIndex);
        structureVoxelCount = structureVoxelCounts(structure);
        
        %- Underdose constraint
        underdoseCoefficient = sqrt(underdoseWeights(structure) / ...
                                    structureVoxelCount);
        underdoseConstraintIndex = v;
        underdoseVarIndex = underdoseConstraintIndex;
        rows(entryIndex) = underdoseConstraintIndex;
        cols(entryIndex) = underdoseVarIndex;
        vals(entryIndex) = -1;
        entryIndex = entryIndex + 1;
        b(underdoseConstraintIndex) = -minDoses(voxelIndex) * ...
                                      underdoseCoefficient;
        
        %- Overdose constraint
        overdoseCoefficient = sqrt(overdoseWeights(structure) / ...
                                   structureVoxelCount);
        overdoseConstraintIndex = underdoseConstraintIndex + ...
                                  numActiveVoxels;
        overdoseVarIndex = overdoseConstraintIndex;
        rows(entryIndex) = overdoseConstraintIndex;
        cols(entryIndex) = overdoseVarIndex;
        vals(entryIndex) = -1;
        entryIndex = entryIndex + 1;
        b(overdoseConstraintIndex) = maxDoses(voxelIndex) * ...
                                     overdoseCoefficient;

        %- Add dosage variables, which are weighted sums of durations
        for shot = 1:numShots
            durationVarIndex = 2*numActiveVoxels + shot;
            durationVarCoefficient = doseRates{shot}(voxelIndex);

            %- Underdose constraint
            rows(entryIndex) = underdoseConstraintIndex;
            cols(entryIndex) = durationVarIndex;
            vals(entryIndex) = -durationVarCoefficient * ...
                               underdoseCoefficient;
            entryIndex = entryIndex + 1;
            
            %- Overdose constraint
            rows(entryIndex) = overdoseConstraintIndex;
            cols(entryIndex) = durationVarIndex;
            vals(entryIndex) = durationVarCoefficient * ...
                               overdoseCoefficient;
            entryIndex = entryIndex + 1;
        end
    end
    A = sparse(rows, cols, vals, numActiveVoxels*2, ...
               numActiveVoxels*2 + numShots);
    
    %- Run quadratic program solver
    fprintf('Starting quadprog.\n');
    sol = quadprog(H, f, A, b, [], [], f, []);
    
    %- Extract durations from solution
    durations = zeros(numIsocenters, numCollimatorSizes);
    for shot = 1:numShots
        durationVarIndex = numActiveVoxels*2 + shot;
        durations(shot) = sol(durationVarIndex);
    end
    
    plan.optimalDurations = durations;
    plan.doseRates = doseRates;
    plan.minDoses = minDoses;
    plan.maxDoses = maxDoses;
    
    plan.totalDoses = zeros(gridSize);
    for c = 1:numIsocenters
        for s = 1:numCollimatorSizes
            duration = plan.optimalDurations(c, s);
            plan.totalDoses = plan.totalDoses + doseRates{c, s} * duration;
        end
    end
end

