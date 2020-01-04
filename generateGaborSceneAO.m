function [scene, I, linearizedI] = generateGaborSceneAO(display, coltype, ort, sf, cont)

% coltype: 1 -- color, 0 -- black white
% ort: 1 -- vertical, 0 -- horizontal

%% HERE'S THE ACTUAL MAKE-GRAITNG PART
TCA = [0, 0]; 

load('gamma_data.mat'); 
% load('LuminousEfficiencyFunction.mat')
inverted_gamma_params = GammaCorrection(gamma_data);
% % RGLuminanceScaling = ScaleRGIntensity(Wavelength, LuminousEfficiencyFunction, R0Intensity, G0Intensity, R255Intensity, G255Intensity);
RGLuminanceScaling = 1;
% 
% BackgroundCol = [0.5, 0.5, 0];
% correctedBackgroundColBit = CorrectGammaBitRG(BackgroundCol, inverted_gamma_params, RGLuminanceScaling);

FOV = 1; %1.5 * 0.5; 
StimulusSizePerDegree = 800; %553;
StimulusSizeInPix = StimulusSizePerDegree * FOV;

GaussianConstant = 4;
current_SpatialFrequency_pixel = sf / StimulusSizePerDegree;

stim = makeGrating (StimulusSizeInPix, ort * 90, current_SpatialFrequency_pixel, 0, GaussianConstant, TCA, coltype);

if length(cont) == 1
    stim_plane_red = cont * stim.R / 2 + 0.5;
    stim_plane_green = cont * stim.G / 2 + 0.5;
elseif length(cont) == 2
    stim_plane_red = cont(1) * stim.R / 2 + 0.5;
    stim_plane_green = cont(2) * stim.G / 2 + 0.5;
else 
    error('Too many values for the image contrast!');
end 
stim = repmat(zeros(size(stim_plane_red)), [1 1 3]);
stim(:,:,1) = stim_plane_red;
stim(:,:,2) = stim_plane_green;
stim = CorrectGammaBitTXTRG(stim, inverted_gamma_params, RGLuminanceScaling);

[scene, I, linearizedI] = sceneFromFile(stim, 'rgb', [], display);

% figure; 
% imshow(stim);

end
