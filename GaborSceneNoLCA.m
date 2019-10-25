clear all; 
close all;

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

theOI_control = oiCreate('wvf human');


%% Create the scene
% presentationDisplay = displayCreate('AOSim-Seattle');
presentationDisplay = displayCreate('AOSim-Seattle_Gabor');
% ppd = 553; using 'AOSim-Seattle_Gabor.mat' - altered d.dpi such that
% (d.dpi/2.54)* (1/((2*atan(1/(2*d.dist*100)))*(180/pi))) == 553 
% from 'AOSim-Seattle.mat'

scene_sample = generateGaborSceneAO(presentationDisplay, 1, 1, 1, 1); % just to get the fov for mosaic generation 
sceneFov = 0.77;%sceneGet(scene_sample, 'fov');
% ok, if this value is smaller than that size of the scene, the mosaic
% generation gets error. Here as a quick remedy artificially giving a
% slightly higer fov. 


%% Generate a hexagonal cone mosaic with ecc-based cone quantal efficiency
KLMSdensity = {[0.0 0.5 0.5 0.0]', [0.0 0.8 0.2 0.0]', [0.0 0.2 0.8 0.0]'}; 

for mos = 1:length(KLMSdensity)
    
    theMosaic = coneMosaicHex(5, ...               % hex lattice sampling factor
       'fovDegs', sceneFov, ...                    % match mosaic width to stimulus size 
       'eccBasedConeDensity', true, ...            % cone density varies with eccentricity
       'eccBasedConeQuantalEfficiency', true, ...  % cone quantal efficiency varies with eccentricity
       'integrationTime', 10/1000, ...             % 30 msec integration time
       'maxGridAdjustmentIterations', 50, ...
       'spatialDensity', KLMSdensity{mos});        % terminate iterative lattice adjustment after 50 iterations
   
    savename = ['mosaic_cond', num2str(mos), '.mat']; 
    save(savename, 'theMosaic'); 
   
    nSF = 7; 
    nContrast = 6; 
    coltype_set = {[0 0], [1 1], [0 1]};     % isochromatic, isoluminant, isochromatic vs. isoluminant
    ort_set = {[0 1], [0 1], [1 1]};         % diff orts, diff orts, same ort
    sf_set = linspace(5, 50, nSF);             % for each of the three above
    contrast_set = linspace(0.001, 0.05, nContrast); % for each of the three above

    condIdPerColtype = [floor((0:nSF*nContrast-1)/nContrast)' + 1, mod(0:nSF*nContrast-1, nContrast)' + 1]; 
    
    for oi = 1:2
        if oi == 2
            theOI = theOI_control; 
            disp('Now using the control OI.'); 
        end 
        
        for exp = 1:length(coltype_set)
            this_coltype = coltype_set{exp}; 
            this_ort     = ort_set{exp}; 
            fprintf('Experimental condition %d. \n', exp); 
            
            result = []; 
            
            for cnd = 1:size(condIdPerColtype, 1)
                close all; 
                
                this_sf       = sf_set(condIdPerColtype(cnd, 1)); 
                this_contrast = contrast_set(condIdPerColtype(cnd, 2)); 
                
                %% ------ LOOP FOR EACH CONDITION FROM HERE ------ %%
                scene1 = generateGaborSceneAO(presentationDisplay, this_coltype(1), this_ort(1), this_sf, this_contrast); % inputs: (display, coltype, ort, sf, contrast)
                scene2 = generateGaborSceneAO(presentationDisplay, this_coltype(2), this_ort(2), this_sf, this_contrast);
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
                
                
                %% Compute some instances of cone mosaic excitations
                nInstancesNum = 64;
                % Zero fixational eye movements
                emPath = zeros(nInstancesNum, 1, 2);
                % Compute mosaic excitation responses
                coneExcitationsCond1 = theMosaic.compute(theOIscene1, 'emPath', emPath);
                coneExcitationsCond2 = theMosaic.compute(theOIscene2, 'emPath', emPath);
                
                % % Compute the mean response across all instances
                % meanConeExcitation = mean(coneExcitationCond1,1);
                % visualizeConeMosaicResponses(theMosaic, coneExcitationCond1, 'R*/cone/tau');
                
                percentCorrect = svm_pca(theMosaic, coneExcitationsCond1, coneExcitationsCond2);
                
                svm_result(cnd) = percentCorrect; 
                fprintf('SF %f, Contrast %f: %f \n', [this_sf, this_contrast, percentCorrect]); 
                
            end
            result{oi} = [result; reshape(svm_result, nSF, nContrast)];
        end
        Result{oi} = result; 
    end
    save(savename, 'Result', '-a'); 
end