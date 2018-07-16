clear; clc; close all;
%% Set up initial parameters
N = 50; R = N/2;
gridSize = [N N]; % = [N N N];
coordinates = voxelgrid.CoordinateData(gridSize);

brainCenter = floor(gridSize / 2) + 1;
brainRadius = R * 0.95;
tumorCorner1 = brainCenter - gridSize/4;
tumorCorner2 = brainCenter + gridSize/4;

numIsocenters = 20;
collimatorSizes = R * [0.05 0.1 0.25];
alpha = 0.8; % closeness acceptance factor
f = 0.15; % relaxation factor

% Structure 1 - Healthy Tissue, Structure 2 - Tumor
healthyTissueMinDose = 0;
healthyTissueMaxDose = 1;
healthyTissueUnderdosingWeight = 0;
healthyTissueOverdosingWeight = 1;
tumorMinDose = 2;
tumorMaxDose = 5;
tumorUnderdosingWeight = 1;
tumorOverdosingWeight = 1; 
%% Initialize structure assignments
healthyTissue = zeros(gridSize);
healthyTissue = voxelgrid.setSphere(healthyTissue, coordinates, ...
                                    brainCenter, brainRadius, 1);
tumor = zeros(gridSize);
tumor = voxelgrid.setAABB(tumor, coordinates, ...
                          tumorCorner1, tumorCorner2, 1);
healthyTissue = healthyTissue - tumor;
structureAssignments = 1*healthyTissue + 2*tumor;
%% Display structure assignments
voxelgrid.displayGrid(structureAssignments);
%% Generate isocenters
isocenters = planning.genIsocenters(tumor, coordinates, numIsocenters, ...
                                    collimatorSizes, alpha, f);
%% Display isocenters (best in 2D)
temp = structureAssignments;
numStructures = max(temp(:));
isocenterLabel = numStructures + 1;
for isocenter = 1:size(isocenters,1)
    loc = num2cell(isocenters(isocenter,:));
    temp(loc{:}) = isocenterLabel;
end
voxelgrid.displayGrid(temp);
%% Calculate dose rates
doseRates = cell(length(isocenters), length(collimatorSizes));
for i = 1:length(isocenters)
    for j = 1:length(collimatorSizes)
        c = isocenters(i,:);
        s = collimatorSizes(j);
        doseRates{i,j} = planning.getDoseRate(coordinates, c, s);
    end
end
%% Generate treatment plan
minDosePerStructure = [healthyTissueMinDose tumorMinDose];
maxDosePerStructure = [healthyTissueMaxDose tumorMaxDose];
underdosingWeights = [healthyTissueUnderdosingWeight, ...
                      tumorUnderdosingWeight];
overdosingWeights = [healthyTissueOverdosingWeight, ...
                     tumorOverdosingWeight]; 

plan = planning.genTreatmentPlan(structureAssignments, ...
                                 minDosePerStructure, ...
                                 maxDosePerStructure, ...
                                 underdosingWeights, overdosingWeights, ...
                                 doseRates);
%% Display treatment plan
voxelgrid.displayGrid(plan.totalDoses);
%% Display volume for which no penalties were applied (excluding unassigned voxels):
satisfied = and(plan.totalDoses < plan.maxDoses, ...
                plan.totalDoses > plan.minDoses);
voxelgrid.displayGrid(satisfied);