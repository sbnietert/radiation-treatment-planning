clear; clc; close all;
%% Set up initial parameters
gridSizeMax = 100;
R = gridSizeMax/2;
numIsocenters = 20;
collimatorSizes = R/5 * [0.05 0.1 0.25];
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
minDosePerStructure = [healthyTissueMinDose tumorMinDose];
maxDosePerStructure = [healthyTissueMaxDose tumorMaxDose];
underdosingWeights = [healthyTissueUnderdosingWeight, ...
                      tumorUnderdosingWeight];
overdosingWeights = [healthyTissueOverdosingWeight, ...
                     tumorOverdosingWeight];
%% Draw structures
brainScan = imread('data/brain tumor.jpg');
imshow(brainScan);
title('Select ROIs - Type Enter into Console to Continue');

brainPoints = load('data/brainPoints.mat');
brainPoints = brainPoints.brainPoints;
brainPoly = impoly(gca, brainPoints);
setColor(brainPoly, 'blue');

tumorPoints = load('data/tumorPoints.mat');
tumorPoints = tumorPoints.tumorPoints;
tumorPoly = impoly(gca, tumorPoints);
setColor(tumorPoly, 'red');
pause
brainSurface = fliplr(brainPoly.getPosition);
tumorSurface = fliplr(tumorPoly.getPosition);
close all;
%% Create voxel grid
imageSize = size(brainScan);
gridSize = imageSize(1:2);
scaleFactor = gridSizeMax / max(gridSize);
gridSize = floor(gridSize * scaleFactor);
coordinates = voxelgrid.CoordinateData(gridSize);
brainSurface = brainSurface * scaleFactor;
tumorSurface = tumorSurface * scaleFactor;
%% Initialize structure assignments
healthyTissue = zeros(gridSize);
healthyTissue = voxelgrid.setSurfaceEnclosedVolume(healthyTissue, ...
                                                   coordinates, ...
                                                   brainSurface, 1);
tumor = zeros(gridSize);
tumor = voxelgrid.setSurfaceEnclosedVolume(tumor, coordinates, ...
                                           tumorSurface, 1);
healthyTissue = healthyTissue - tumor;

structureAssignments = 1*healthyTissue + 2*tumor;
%% Display structure assignments
voxelgrid.displayGrid(structureAssignments);
%% Generate isocenters
isocenters = planning.genIsocenters(tumor, coordinates, numIsocenters, ...
                                    collimatorSizes, alpha, f);
%% Display isocenters
temp = structureAssignments;
numStructures = max(temp(:));
isocenterLabel = numStructures + 1;
for isocenter = 1:size(isocenters,1)
    loc = num2cell(isocenters(isocenter,:));
    temp(loc{:}) = isocenterLabel;
end
voxelgrid.displayGrid(temp);
%% Calculate dose rates
doseRates = cell(size(isocenters, 1), length(collimatorSizes));
for i = 1:size(isocenters, 1)
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