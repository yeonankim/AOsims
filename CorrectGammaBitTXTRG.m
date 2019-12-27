function correctedColBit = CorrectGammaBitTXTRG(Col, inverted_gamma_params, RGLuminanceScaling, ColorChannel)

  RParameter = [inverted_gamma_params(1).Parameter 0];
  GParameter = [inverted_gamma_params(2).Parameter 0];
  
  if ~ exist('RGLuminanceScaling')
    RGLuminanceScaling = 1;
  elseif isempty(RGLuminanceScaling)
      RGLuminanceScaling = 1;
  end

  R = Col(:,:,1);
  G = Col(:,:,2) / RGLuminanceScaling;
  B = Col(:,:,3);

  correctedR = InvertGammaExtP(RParameter, 255, R);
  correctedG = InvertGammaExtP(GParameter, 255, G);

  correctedColBit = round(cat(3, correctedR,correctedG,B * 255));

  if exist('ColorChannel')
    if strcmp(ColorChannel, 'R')
      correctedColBit(:,:,2) = 0;
    elseif strcmp(ColorChannel, 'G')
      correctedColBit(:,:,1) = 0;
    end
  end

end
