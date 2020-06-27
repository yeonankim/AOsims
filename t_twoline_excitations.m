% Basic tutorial of the UW "Two line Scene" using ISETBio No-LCA optics
%
% Description:
%   ISETBio's wavefront optics normally includes wavelength-dependent
%   defocus to model longitudinal chromatic aberration (LCA). This tutorial shows
%   how to construct an optical image object that does not do that.  Useful
%   for modeling experiments where LCA has been corrected for.
%
%   Retinal images are produced for diffraction limited seeing, with and
%   without LCA.  For the case with, stimuli at the accommodated wavelength
%   are diffraction limited, while at other wavelengths a defocus terms is
%   added to the optics specification according to typical human LCA.
%
%   You can adjust the accommodated wavelength and see the effect of LCA in
%   one of the computed retinal images and not the other.  The test scene
%   used employs wavelengths of 530 and 660.

% History:
%  05/01/19  dhb  Version with LCA defeated.
%  06/27/19  drc  Move make_optics() to another file, here just scene.

[theOIWithLca,theOINoLca] = make_optics();

presentationDisplay = displayCreate('AOSim-Seattle_SPDcorrected_Scaled');

% Generate TwoLine Scene
scene = generateTwoLineScene(presentationDisplay,1,20); % 1=RG
sceneFov = sceneGet(scene, 'fov');

%% Compute and visualize the retinal images with and without LCA
theOINoLca = oiCompute(theOINoLca, scene);
visualizeOpticalImage(theOINoLca, 'displayRadianceMaps', false, ...
    'displayRetinalContrastProfiles', false);

theOIWithLca = oiCompute(theOIWithLca, scene);
% visualizeOpticalImage(theOIWithLca, 'displayRadianceMaps', false, ...
%   'displayRetinalContrastProfiles', false);

this_KLMSdensity=[0.0 10/16 5/16 1/16]';              % for example: 2L 1M 0S
theMosaic = coneMosaicHex(5, ...               % hex lattice sampling factor
    'fovDegs', sceneFov, ...                    % match mosaic width to stimulus size
    'eccBasedConeDensity', true, ...            % cone density varies with eccentricity
    'eccBasedConeQuantalEfficiency', true, ...  % cone quantal efficiency varies with eccentricity
    'integrationTime', 10/1000, ...             % 30 msec integration time
    'maxGridAdjustmentIterations', 50, ...
    'spatialDensity', this_KLMSdensity);        % terminate iterative lattice adjustment after 50 iterations

%% Compute some instances of cone mosaic excitations
nInstancesNum = 10;
% Zero fixational eye movements
emPath = zeros(nInstancesNum, 1, 2);
% Compute mosaic excitation responses
coneExcitations = theMosaic.compute(theOINoLca, 'emPath', emPath);
% Compute the mean response across all instances
meanConeExcitation = mean(coneExcitations,1);

%coneExcitations = theMosaic.compute(theOINoLca);

visualizeConeMosaicResponses(theMosaic, coneExcitations, 'R*/cone/tau');
%%This causes an error for DRC on the second column of plot 6/27/2020


