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
presentationDisplay = displayCreate('AOSim-Seattle');

scene_sample = generateGaborSceneAO(presentationDisplay, 1, 1, 1, 1); % just to get the fov for mosaic generation 
sceneFov = 1.05;%sceneGet(scene_sample, 'fov');
% ok, if this value is smaller than that size of the scene, the mosaic
% generation gets error. Here as a quick remedy artificially giving a
% slightly higer fov. 


%% Generate a hexagonal cone mosaic with ecc-based cone quantal efficiency
KLMSdensity = {[0.0 0.5 0.5 0.0]', [0.0 0.9 0.1 0.0]'};%, [0.0 0.2 0.8 0.0]'};
nSF = 3;
nContrast = 6;
coltype_set = {[0 0], [1 1], [0 1]};     % isochromatic, isoluminant, isochromatic vs. isoluminant
ort_set = {[0 1], [0 1], [1 1]};         % diff orts, diff orts, same ort
sf_set = [10 30 60]; %linspace(5, 50, nSF);             % for each of the three above
contrast_set = 10.^linspace(log10(0.001), log10(0.03), nContrast); %linspace(0.001, 0.05, nContrast); % for each of the three above

% Making dir to save all the simulation results
resultdir = 'Results';
if ~isfolder(resultdir)
    mkdir(resultdir);
end

for mos = 1:length(KLMSdensity)
    
    this_KLMSdensity = KLMSdensity{mos};
    
    % Making dir to save cone excitation instances
    mosaicdir = fullfile(resultdir, ['mosaicCond', num2str(mos)]);
    if ~isfolder(mosaicdir)
        mkdir(mosaicdir);
    end

    savename_mosaic = fullfile(mosaicdir, ['mosaic_L', num2str(this_KLMSdensity(2)*10), ...,
                                                  'M', num2str(this_KLMSdensity(3)*10), ...,
                                                  'S', num2str(this_KLMSdensity(4)*10), '.mat']);
                                              
    if isfile(savename_mosaic)
        load(savename_mosaic);
        fprintf('Loading file: %s \n', savename_mosaic);
    else
        theMosaic = coneMosaicHex(5, ...               % hex lattice sampling factor
            'fovDegs', sceneFov, ...                    % match mosaic width to stimulus size
            'eccBasedConeDensity', true, ...            % cone density varies with eccentricity
            'eccBasedConeQuantalEfficiency', true, ...  % cone quantal efficiency varies with eccentricity
            'integrationTime', 10/1000, ...             % 30 msec integration time
            'maxGridAdjustmentIterations', 50, ...
            'spatialDensity', this_KLMSdensity, ...     % terminate iterative lattice adjustment after 50 iterations
            'noiseFlag', 'none');
        save(savename_mosaic, 'theMosaic', 'this_KLMSdensity');
    end 
    
    
    % Making dir to save cone excitation instances
    conerespdir = fullfile(mosaicdir, 'ConeExitationInstances');
    if ~isfolder(conerespdir)
        mkdir(conerespdir);
    end

    condIdPerColtype = [floor((0:nSF*nContrast-1)/nContrast)' + 1, mod(0:nSF*nContrast-1, nContrast)' + 1]; 
    
    for oi = 1 %1:2 **Running only for AO setup
        if oi == 2
            theOI = theOI_control; 
            disp('Now using the control OI.'); 
        end 
        
        for exp = 1:length(coltype_set) %**Running only for experiment 3
            
            this_coltype = coltype_set{exp}; 
            this_ort     = ort_set{exp}; 
            fprintf('Experimental condition %d. \n', exp); 
            
            parfor cnd = 1:size(condIdPerColtype, 1)
                close all; 
                
                this_sf       = sf_set(condIdPerColtype(cnd, 1)); 
                this_contrast = contrast_set(condIdPerColtype(cnd, 2)); 
                
                savename_coneresp = fullfile(conerespdir, ['coneExcitation_noiseOff_exp', num2str(exp), '_SF_', num2str(this_sf), '_contr_', num2str(this_contrast), '.mat']); 
                
                if isfile(savename_coneresp)
                    fprintf('This cone excitation instance already exists. Skipping... \n');
                else
                    
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
                    %                 nInstancesNum = 1024;
                    %                 % Zero fixational eye movements
                    %                 emPath = zeros(nInstancesNum, 1, 2);
                    %                 % Compute mosaic excitation responses
                    %                 coneExcitationsCond1 = theMosaic.compute(theOIscene1, 'emPath', emPath);
                    %                 coneExcitationsCond2 = theMosaic.compute(theOIscene2, 'emPath', emPath);
                    
                    coneExcitationsCond1 = theMosaic.compute(theOIscene1);
                    coneExcitationsCond2 = theMosaic.compute(theOIscene2);
                    
                    saveparfor(savename_coneresp, coneExcitationsCond1, coneExcitationsCond2)
                    
                end 
                
%                 % % Compute the mean response across all instances
%                 % meanConeExcitation = mean(coneExcitationCond1,1);
%                 % visualizeConeMosaicResponses(theMosaic, coneExcitationCond1, 'R*/cone/tau');
%                 
%                 percentCorrect = svm_pca(theMosaic, coneExcitationsCond1, coneExcitationsCond2);
%                 
%                 fprintf('SF %f, Contrast %f: %f \n', [this_sf, this_contrast, percentCorrect]);                 
%                 Result{oi, exp}(mod(cnd-1, nContrast)+1, floor((cnd-1)/nContrast)+1) = percentCorrect; 
%                 
%                 save(savename, 'Result', '-append'); 
%                 
%                 savename_coneInst = [parentdir, '/mosaicCond_', num2str(mos), '_oiCond_', num2str(oi), '_exp_', num2str(exp)]; 
%                 save(savename_coneInst, 'coneExcitationsCond1', 'coneExcitationsCond2', '-v7.3'); 
            end
        end
    end    
end



function saveparfor(savename_coneresp, coneExcitationsCond1, coneExcitationsCond2)

fprintf('Saving %s. \n', savename_coneresp);
save(savename_coneresp, 'coneExcitationsCond1', 'coneExcitationsCond2');

end 