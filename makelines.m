% This function, originally from Ram Sabesan in UW Seattle, creates the
% stimulus bits for AO Vision Simulator "two-Line" experiment.
%
% Stimuli are two side-by-side small vertical lines (arbitrary thickness and
% spacing), which can be one of four color combinations: RG, GR, YY, or
% a single yellow. Yellow is 0.5R + 0.5G.
%
% Returned datatype is a structure with separate images planes for R and G.

function  img = makelines (imageSize, TCA, stimcode, spacing,line_thickness,line_height,background_level)

%   %testl
%   imageSize = 401;
%   stimcode = 1;
%   spacing = 10;
%   line_thickness = 4;
%   line_height    = 200;
%   %test end
%   

  img_center = floor([imageSize/2,imageSize/2]); 
  img        = ones(imageSize)*background_level;
  img_red    = img;
  img_green  = img;
  
  switch(stimcode)
  
      case 1 %RG
          img_red(img_center - line_height/2:img_center + line_height/2,...
              img_center - spacing/2 - line_thickness : img_center - spacing/2-1) = 1; 
          
          img_green(img_center - line_height/2:img_center + line_height/2,...
              img_center + spacing/2 : img_center + spacing/2 + line_thickness-1) = 1; 
         
  
      case 2 %GR
          
          img_green(img_center - line_height/2:img_center + line_height/2,...
              img_center - spacing/2 - line_thickness : img_center - spacing/2-1) = 1; 
          
          img_red(img_center - line_height/2:img_center + line_height/2,...
              img_center + spacing/2 : img_center + spacing/2 + line_thickness-1) = 1; 
                            
      case 3 %YY
          
          img_green(img_center - line_height/2:img_center + line_height/2,...
              img_center - spacing/2 - line_thickness : img_center - spacing/2-1) = 0.5; 
          
          img_green(img_center - line_height/2:img_center + line_height/2,...
              img_center + spacing/2 : img_center + spacing/2 + line_thickness-1) = 0.5;
          
          
          img_red(img_center - line_height/2:img_center + line_height/2,...
              img_center - spacing/2 - line_thickness : img_center - spacing/2-1) = 0.5; 
          
          img_red(img_center - line_height/2:img_center + line_height/2,...
              img_center + spacing/2 : img_center + spacing/2 + line_thickness-1) = 0.5; 
          
          
      case 4 % YY-abut
          spacing = 0;
          
          img_green(img_center - line_height/2:img_center + line_height/2,...
              img_center - spacing/2 - line_thickness : img_center - spacing/2-1) = 0.5; 
          
          img_green(img_center - line_height/2:img_center + line_height/2,...
              img_center + spacing/2 : img_center + spacing/2 + line_thickness-1) = 0.5;
          
          
          img_red(img_center - line_height/2:img_center + line_height/2,...
              img_center - spacing/2 - line_thickness : img_center - spacing/2-1) = 0.5; 
          
          img_red(img_center - line_height/2:img_center + line_height/2,...
              img_center + spacing/2 : img_center + spacing/2 + line_thickness-1) = 0.5; 
  end
   
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
