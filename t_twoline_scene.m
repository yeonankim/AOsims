% Show how to create an oi structure with no LCA
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

%% Compute and visualize the retinal images with and without LCA
theOINoLca = oiCompute(theOINoLca, scene);
visualizeOpticalImage(theOINoLca, 'displayRadianceMaps', false, ...
    'displayRetinalContrastProfiles', false);

theOIWithLca = oiCompute(theOIWithLca, scene);
visualizeOpticalImage(theOIWithLca, 'displayRadianceMaps', false, ...
    'displayRetinalContrastProfiles', false);
