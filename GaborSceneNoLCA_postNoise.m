clear all;
close all;

rng;

[theOI_control,theOI]=make_optics();

%% Create the scene
presentationDisplay = displayCreate('AOSim-Seattle_SPDcorrected_Scaled');

scene_sample = generateGaborSceneAO(presentationDisplay, 1, 1, 1, 1); % just to get the fov for mosaic generation
sceneFov = 1.1;%sceneGet(scene_sample, 'fov');
% ok, if this value is smaller than that size of the scene, the mosaic
% generation gets error. Here as a quick remedy artificially giving a
% slightly higer fov.


%% Generate a hexagonal cone mosaic with ecc-based cone quantal efficiency
%-------------------------------------------------%
resultdir = 'Results';
foldername = 'ConeExcitationInstances_SPDcorrectedScaled_60PCA_1024Instances';
KLMSdensity = {[0 13/14 1/14 0]', [0 5/6 1/6 0]', [0 2.8/3.8 1/3.8 0]', [0 1.8/2.8 1/2.8 0]', [0.0 0.5 0.5 0.0]', [0 1/2.8 1.8/2.8 0]', [0 1/3.8 2.8/3.8, 0]', [0 1/6 5/6 0]', [0 1/14 13/14 0]'};
coltype_set = {[0 0], [1 1], [0 1]};     % isochromatic, isoluminant, isochromatic vs. isoluminant
ort_set = {[0 1], [0 1], [1 1]};         % diff orts, diff orts, same ort
sf_set = [4, 8, 16, 24, 32, 48, 64]; 
nSF = length(sf_set); %[10 30 60]; %linspace(5, 50, nSF);             % for each of the three above
nContrast = 15; contrast_set = 10.^linspace(log10(0.001), log10(0.08), nContrast); %10.^linspace(log10(0.03), log10(0.1), nContrast); %linspace(0.001, 0.05, nContrast); % for each of the three above
nSVMrep = 1; 
%-------------------------------------------------%


% Making dir to save all the simulation results
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
        if ~strcmp(theMosaic.noiseFlag, 'none')
            fprintf('The noiseFlag of this mosaic is set to %s. Turning off the noise. /n', theMosaic.noiseFlag);
            theMosaic.noiseFlag = 'none';
        end
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
    conerespdir = fullfile(mosaicdir, foldername);
    if ~isfolder(conerespdir)
        mkdir(conerespdir);
        mkdir(fullfile(conerespdir, 'noisyInstances'));
    end
    
    condIdPerColtype = [floor((0:nSF*nContrast-1)/nContrast)' + 1, mod(0:nSF*nContrast-1, nContrast)' + 1];
    
    for oi = 1:2
        % For oi==1, do LCA off, everything OFF
        % For oi==2, do LCA on, normal wvf
        if oi == 2
            theOI = theOI_control;
            disp('Now using the control OI.');
        end
        
        for exp = 1:length(coltype_set) 
            
            this_coltype = coltype_set{exp};
            this_ort     = ort_set{exp};
            fprintf('Experimental condition %d. \n', exp);
            
            parfor cnd = 1:size(condIdPerColtype, 1)
                close all;
                
                this_sf       = sf_set(condIdPerColtype(cnd, 1));
                this_contrast = contrast_set(condIdPerColtype(cnd, 2));
                
                savename_coneresp = fullfile(conerespdir, ['coneExcitation_noiseOff_oi', num2str(oi), '_exp', num2str(exp), '_SF_', num2str(this_sf), '_contr_', num2str(this_contrast), '.mat']);
                [path, fname, ext] = fileparts(savename_coneresp);
                savename_coneresp_instances = fullfile(path, 'noisyInstances', ['noisyInstances_', fname, ext]);
                
                if isfile(savename_coneresp)
                    fprintf('This cone excitation instance already exists. Skip computing... \n');
                    [runsvm, nrunadd] = checkforsvm(savename_coneresp, nSVMrep);
                    if runsvm
                        if isfile(savename_coneresp_instances)
                            stemp_instances = load(savename_coneresp_instances);
                            coneExcitationCond1_noisyInstances = stemp_instances.coneExcitationCond1_noisyInstances;
                            coneExcitationCond2_noisyInstances = stemp_instances.coneExcitationCond2_noisyInstances;
                        else
                            stemp = load(savename_coneresp);
                            [coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances] = coneExcitationAddNoise(stemp.coneExcitationsCond1, stemp.coneExcitationsCond2);
                            saveparfor_instances(savename_coneresp_instances, coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances);
                        end
                        SVMpercentCorrect = [];
                        for rep = 1:nrunadd
                            SVMpercentCorrect = [SVMpercentCorrect, svm_pca(theMosaic, coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances)];
                        end
                        saveparfor_svm(savename_coneresp, SVMpercentCorrect);
                    end
                    
                else
                    %% ------ LOOP FOR EACH CONDITION FROM HERE ------ %%
                    scene1 = generateGaborSceneAO(presentationDisplay, this_coltype(1), this_ort(1), this_sf, this_contrast); % inputs: (display, coltype, ort, sf, contrast)
                    scene2 = generateGaborSceneAO(presentationDisplay, this_coltype(2), this_ort(2), this_sf, this_contrast);
%                     scene2 = generateGaborSceneAO(presentationDisplay, this_coltype(2), this_ort(2), this_sf, 0);
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
                    
                    
                    %% Compute the mean cone excitation responses to the stimulus
                    coneExcitationsCond1 = theMosaic.compute(theOIscene1);
                    coneExcitationsCond2 = theMosaic.compute(theOIscene2);
                    saveparfor(savename_coneresp, coneExcitationsCond1, coneExcitationsCond2); 
                    
                    randomSeed = rng;
                    saveparfor_rseed(savename_coneresp, randomSeed); 
                    
                    [coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances] = coneExcitationAddNoise(coneExcitationsCond1, coneExcitationsCond2);
                    saveparfor_instances(savename_coneresp_instances, coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances);
                    
                    SVMpercentCorrect = [];
                    for rep = 1:nSVMrep
                        SVMpercentCorrect = [SVMpercentCorrect, svm_pca(theMosaic, coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances)];
                    end
                    saveparfor_svm(savename_coneresp, SVMpercentCorrect);
                    
                end

            end
        end
    end
end

function [coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances] = coneExcitationAddNoise(coneExcitationsCond1, coneExcitationsCond2)

% Compute some noisy instances of cone mosaic excitations
nInstances = 1024;
[nr1, nc1] = size(coneExcitationsCond1);
[nr2, nc2] = size(coneExcitationsCond2);
coneExcitationCond1_noisyInstances = nan(nInstances, nr1, nc1); 
coneExcitationCond2_noisyInstances = nan(nInstances, nr2, nc2); 
for i = 1:nInstances
    coneExcitationCond1_noisyInstances(i,:,:) = poissrnd(coneExcitationsCond1);
    coneExcitationCond2_noisyInstances(i,:,:) = poissrnd(coneExcitationsCond2);
end

end


function saveparfor(savename_coneresp, coneExcitationsCond1, coneExcitationsCond2)

fprintf('Saving %s. \n', savename_coneresp);
save(savename_coneresp, 'coneExcitationsCond1', 'coneExcitationsCond2');

end

function saveparfor_instances(savename_coneresp_instances, coneExcitationCond1_noisyInstances, coneExcitationCond2_noisyInstances)

save(savename_coneresp_instances, 'coneExcitationCond1_noisyInstances', 'coneExcitationCond2_noisyInstances', '-v7.3');

end

function saveparfor_rseed(savename_coneresp, randomSeed)

save(savename_coneresp, 'randomSeed', '-append');

end

function saveparfor_svm(savename_coneresp, SVMpercentCorrect)

fprintf('%s: %.2f \n', savename_coneresp, mean(SVMpercentCorrect));
if ismember(who('-file', savename_coneresp), 'SVMpercentCorrect')
    SVMpercentCorrect_prev = load(savename_coneresp, 'SVMpercentCorrect');
    SVMpercentCorrect = [SVMpercentCorrect_prev, SVMpercentCorrect];
end
save(savename_coneresp, 'SVMpercentCorrect', '-append');

end

function [runsvm, nrunadd] = checkforsvm(savename_coneresp, nSVMrep)
if ismember(who('-file', savename_coneresp), 'SVMpercentCorrect')
    SVMpercentCorrect = load(savename_coneresp, 'SVMpercentCorrect'); 
end 
if length(SVMpercentCorrect) < nSVMrep
    runsvm = true;
    nrunadd = nSVM-length(SVMpercentCorrect);
    if ~ismember(who('-file', savename_coneresp), 'randomSeed')
        randomSeed = rng;
        saveparfor_rseed(savename_coneresp, randomSeed);
    end
    fprintf('%d SVM run exists. Adding %d more runs... \n', length(SVMpercentCorrect), nrunadd);
else
    runsvm = false;
    nrunadd = 0;
    fprintf('%d SVM run exists - this simulation is complete. \n', length(SVMpercentCorrect)); 
end 
end 
