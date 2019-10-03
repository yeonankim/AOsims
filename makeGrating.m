function  img = makeGrating (imageSize, angle, sf, phase, gaussianConstant, TCA, counter_phase_RG)


  if nargin < 4 | isempty(phase)
    phase = 0;
  end
  if nargin<3 | isempty(sf)
    sf = 0.02;
  end
  if nargin<2 | isempty(angle)
    angle = 0;
  end
  if nargin<1 | isempty(imageSize)
    imageSize = 100;
  end

  if ~ exist('counter_phase_RG')
    img = bwGrating(imageSize, angle, sf, phase, gaussianConstant, TCA);
  else
    if counter_phase_RG == 0
      img = bwGrating(imageSize, angle, sf, phase, gaussianConstant, TCA);
    elseif counter_phase_RG == 1
      img = interleavedGrating(imageSize, angle, sf, phase, gaussianConstant, TCA);
    end
  end

end

function img = simpleGrating(imageSize, angle, sf, phase, gaussianConstant)

  x = 1 : imageSize;
  x = x - mean(x);
  [x, y]= meshgrid(x, x);
  angle = angle * pi / 180;
  phase = phase * pi / 180;

  % x = x./max(x(:));
  % y=y./max(y(:));
  % img = (x.^2 + y.^2) <=1; % for alignement

  img = sin(phase + 2 * pi * sf * (cos(angle) * y + sin(angle) * x));
  if exist('gaussianConstant', 'var') & gaussianConstant > 1
    spaceConstant = imageSize / gaussianConstant;
    img = img .* exp(- ((x / spaceConstant) .^ 2)-(( y / spaceConstant) .^ 2));
  end

end

function img = bwGrating(imageSize, angle, sf, phase, gaussianConstant, TCA)

  img = simpleGrating(imageSize, angle, sf, phase, gaussianConstant);
  img_red = img;
  img_green = img;

  img = CorrectTCA(img_red,img_green, TCA);

end

function img = interleavedGrating(imageSize, angle, sf, phase, gaussianConstant, TCA)

  img = simpleGrating(imageSize, angle, sf, phase, gaussianConstant);
  img_red = img;
  img_green = - img;

  img = CorrectTCA(img_red,img_green, TCA);

end

function correctedImage = CorrectTCA(RedImage,GreenImage, TCA)

  horizontlTCA = TCA(1);
  verticalTCA = TCA(2);
  absolute_hozitontalTCA = abs(horizontlTCA);
  absolute_verticalTCA = abs(verticalTCA);

  vertical_crop_pixels = ceil(absolute_verticalTCA / 2);
  horizontal_crop_pixels = ceil(absolute_hozitontalTCA / 2);

  croppedRedImage = RedImage(1 + vertical_crop_pixels : end - vertical_crop_pixels, 1 + horizontal_crop_pixels : end - horizontal_crop_pixels);
  croppedGreenImage = GreenImage(1 + vertical_crop_pixels : end - vertical_crop_pixels, 1 + horizontal_crop_pixels : end - horizontal_crop_pixels);

  if horizontlTCA > 0

    correctedRedImage = [zeros(size(croppedRedImage, 1), absolute_hozitontalTCA) croppedRedImage];
    correctedGreenImage = [croppedGreenImage zeros(size(croppedGreenImage, 1), absolute_hozitontalTCA)];

  elseif horizontlTCA < 0

    correctedGreenImage = [zeros(size(croppedRedImage, 1), absolute_hozitontalTCA) croppedGreenImage];
    correctedRedImage = [croppedRedImage zeros(size(croppedRedImage, 1), absolute_hozitontalTCA)];

  else

    correctedRedImage = croppedRedImage;
    correctedGreenImage = croppedGreenImage;

  end

  if verticalTCA > 0

    correctedRedImage = [zeros(absolute_verticalTCA, size(correctedRedImage, 2)); correctedRedImage];
    correctedGreenImage = [correctedGreenImage; zeros(absolute_verticalTCA, size(correctedGreenImage, 2))];

  elseif verticalTCA < 0

    correctedGreenImage = [zeros(absolute_verticalTCA, size(correctedGreenImage, 2)); correctedGreenImage];
    correctedRedImage = [correctedRedImage; zeros(absolute_verticalTCA, size(correctedRedImage, 2))];

  end

  correctedImage.R = correctedRedImage;
  correctedImage.G = correctedGreenImage;

end
