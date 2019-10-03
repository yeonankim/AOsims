function scene = generateGaborSceneAO(display)
    
% %% Copy from two line example
% % stimcode: 1=RG 2=GR 3=YY 4=Y
% background_level=0.1; % luminance to put in each gun as background
% 
% % All these are in pixels:
% imsize=140;
% %separation=10;
% width=4;
% height=20;
% tca=[0 0]; % Mostly for real experiment, which required subjective adjustment for TCA
% 
% % "Makelines" makes the two stimulus planes
% bits=makelines(imsize,tca,stimcode,separation,width,height,background_level);
% img=zeros(imsize,imsize,3);
% img(:, :, 1) = bits.R;
% img(:, :, 2) = bits.G;


%% HERE'S THE ACTUAL MAKE-GRAITNG PART
TCA = [0, 0]; 

% load('gamma_data.mat')
% load('LuminousEfficiencyFunction.mat')
% inverted_gamma_params = GammaCorrection(gamma_data);
% % RGLuminanceScaling = ScaleRGIntensity(Wavelength, LuminousEfficiencyFunction, R0Intensity, G0Intensity, R255Intensity, G255Intensity);
% RGLuminanceScaling = 1;
% 
% BackgroundCol = [0.5, 0.5, 0];
% correctedBackgroundColBit = CorrectGammaBitRG(BackgroundCol, inverted_gamma_params, RGLuminanceScaling);


FOV = 1.5;

StimulusSizePerDegree = 553*0.5;
StimulusSizeInPix = StimulusSizePerDegree * FOV;
StimulusSizeInDeg = FOV;

% SpatialFrequency = [1 2 4 8 16 24 32 40 48 64 80 96 128 160 192 256]; % cycle / degree
% SpatialFrequency = [1 2 4 8 16 24 32];
SpatialFrequency = [2,7 ];
% StimulusCycles = SpatialFrequency * StimulusSizeInDeg;

GaussianConstant = 4;
current_SpatialFrequency = SpatialFrequency(2);
current_SpatialFrequency_pixel = current_SpatialFrequency * StimulusSizeInDeg / StimulusSizeInPix;

StimulusSequence = [0 1];
% StimulusSequence = StimulusSequence(randperm(length(StimulusSequence)));
% 1 -- color, 0 -- black white

%     orientation = StimulusSequence(randperm(length(StimulusSequence)));
% 1 -- vertical, 0 -- horizontal
orientation = StimulusSequence(randperm(length(StimulusSequence))) * 0;
current_Contrast = 0.5; 

stim1 = makeGrating (StimulusSizeInPix, orientation(1) * 90 + round(rand) * 180, current_SpatialFrequency_pixel, 0, GaussianConstant, TCA, StimulusSequence(1));
stim2 = makeGrating (StimulusSizeInPix, orientation(2) * 90 + round(rand) * 180, current_SpatialFrequency_pixel, 0, GaussianConstant, TCA, StimulusSequence(2));

stim1_plane_red = current_Contrast * stim1.R / 2 + 0.5;
stim1_plane_green = current_Contrast * stim1.G / 2 + 0.5;
stim1 = repmat(zeros(size(stim1_plane_red)), [1 1 3]);
stim1(:,:,1) = stim1_plane_red;
stim1(:,:,2) = stim1_plane_green;
% stim1 = CorrectGammaBitTXTRG(stim1, inverted_gamma_params, RGLuminanceScaling);

stim2_plane_red = current_Contrast * stim2.R / 2 + 0.5;
stim2_plane_green = current_Contrast * stim2.G / 2 + 0.5;
stim2 = repmat(zeros(size(stim2_plane_red)), [1 1 3]);
stim2(:,:,1) = stim2_plane_red;
stim2(:,:,2) = stim2_plane_green;
% stim2 = CorrectGammaBitTXTRG(stim2, inverted_gamma_params, RGLuminanceScaling);


figure; 
subplot(1,2,1); imshow(stim1);
subplot(1,2,2); imshow(stim2);

scene = sceneFromFile(stim1, 'rgb', [], display);
    
end
