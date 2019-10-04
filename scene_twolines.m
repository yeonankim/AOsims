function GenerateTwoLineScene()
stimcode=1; % 1=RG 2=GR 3=YY 4=Y
background_level=0.1; % luminance to put in each gun as background

% All these are in pixels:
imsize=401;
separation=10;
width=4;
height=100;
tca=[0 0]; % Mostly for real experiment, which required subjective adjustment for TCA

% "Makelines" makes the two stimulus planes
bits=makelines(imsize,tca,stimcode,separation,width,height,background_level);
img=zeros(imsize,imsize,3);
img(:, :, 1) = bits.R*10;
img(:, :, 2) = bits.G;

scene = sceneFromFile(img, 'rgb', [], d);