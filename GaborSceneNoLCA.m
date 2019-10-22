clear; close all;

%% Parameters
%
% Make a zero vector of Zernike coefficients to
% represent a diffraction limited pupil function,
% and a few other things.
pupilDiameterMm = 6;
wave = (400:10:700)';
accommodatedWavelength = 530;
zCoeffs = zeros(66,1);

LCAoff = true; 
opticsName = 'human-wvf-withlca'; 


%% Set up wavefront optics object
% Compute pupil function using 'no lca' key/value pair to turn off LCA.
% You can turn it back on to compare the effect.
wvfP = wvfCreate('calc wavelengths', wave, ...
                 'zcoeffs', zCoeffs, ...
                 'measured pupil size', pupilDiameterMm, ...
                 'calc pupil size', pupilDiameterMm, ...
                 'measured wavelength', accommodatedWavelength, ...
                 'name', sprintf('human-%d', pupilDiameterMm));
% Deal with best focus by specifying that the wavefront parameters
% were measured at the wavelength we want to say is in focus. This
% is a little bit of a hack but seems OK for the diffraction limited case
% we're using here.


%% Make optical image object using wvfP and no LCA calc
% Same as above but don't defeat LCA calc
wvfP = wvfComputePupilFunction(wvfP,false,'no lca', LCAoff);
wvfP = wvfComputePSF(wvfP);
theOI = wvf2oi(wvfP);
optics = oiGet(theOI, 'optics');
if LCAoff; opticsName = 'human-wvf-nolca'; end  
optics = opticsSet(optics, 'model', 'shift invariant', 'name', opticsName);
theOI = oiSet(theOI,'optics',optics);

% theOI = oiCreate('wvf human');


%% Create the scene
presentationDisplay = displayCreate('AOSim-Seattle');



%% ------ LOOP FOR EACH CONDITION FROM HERE ------ %%
scene1 = generateGaborSceneAO(presentationDisplay, 0, 50, 0, .05); % inputs: (display, coltype, sf, ort, contrast)
scene2 = generateGaborSceneAO(presentationDisplay, 0, 50, 1, .05); % inputs: (display, coltype, sf, ort, contrast)

% scene = generateTwoLineScene(presentationDisplay,1,20);
% visualizeScene(scene, 'displayContrastProfiles', true); 


%% Compute and visualize the retinal images with and without LCA
theOIscene1 = oiCompute(theOI, scene1);
theOIscene2 = oiCompute(theOI, scene2);

% % Visualize the PSFs and OTFs
% % Visualize the PSF/OTF at 530 (in-focus)
% visualizedSpatialSupportArcMin = 6;
% visualizedSpatialSfrequencyCPD = 120;
% visualizeOptics(theOIscene1, accommodatedWavelength, visualizedSpatialSupportArcMin, visualizedSpatialSfrequencyCPD);
% visualizeOptics(theOIscene2, accommodatedWavelength, visualizedSpatialSupportArcMin, visualizedSpatialSfrequencyCPD);

% % Visualize the optical image as an RGB image and a few spectral slices
% visualizeOpticalImage(theOIscene1, 'displayRadianceMaps', false, ...
%     'displayRetinalContrastProfiles', false);
% visualizeOpticalImage(theOIscene2, 'displayRadianceMaps', false, ...
%     'displayRetinalContrastProfiles', false);


%% Generate a hexagonal cone mosaic with ecc-based cone quantal efficiency
theMosaic = coneMosaicHex(5, ...               % hex lattice sampling factor
   'fovDegs', sceneGet(scene1, 'fov'), ...      % match mosaic width to stimulus size
   'eccBasedConeDensity', true, ...            % cone density varies with eccentricity
   'eccBasedConeQuantalEfficiency', true, ...  % cone quantal efficiency varies with eccentricity
   'integrationTime', 30/1000, ...             % 30 msec integration time
   'maxGridAdjustmentIterations', 50, ...
   'spatialDensity', [0 0.9 0.1 0.0]');        % terminate iterative lattice adjustment after 50 iterations

%% Compute some instances of cone mosaic excitations
nInstancesNum = 32;
% Zero fixational eye movements
emPath = zeros(nInstancesNum, 1, 2);
% Compute mosaic excitation responses
coneExcitationsCond1 = theMosaic.compute(theOIscene1, 'emPath', emPath);
coneExcitationsCond2 = theMosaic.compute(theOIscene2, 'emPath', emPath);

% % Compute the mean response across all instances
% meanConeExcitation = mean(coneExcitationCond1,1);
% visualizeConeMosaicResponses(theMosaic, coneExcitationCond1, 'R*/cone/tau');

svm_pca(theMosaic, coneExcitationsCond1, coneExcitationsCond2)