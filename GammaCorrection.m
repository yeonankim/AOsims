function [inverted_gamma_params, inverted_gamma] = GammaCorrection(gamma_data)

  RBit = [gamma_data(strcmp({gamma_data(:).channel}, 'R')).bit];
  RIntensity = [gamma_data(strcmp({gamma_data(:).channel}, 'R')).intensity];
  [RBit, normedRIntensity] = NormalizeIntensity(RBit, RIntensity);

  GBit = [gamma_data(strcmp({gamma_data(:).channel}, 'G')).bit];
  GIntensity = [gamma_data(strcmp({gamma_data(:).channel}, 'G')).intensity];
  [GBit, normedGIntensity] = NormalizeIntensity(GBit, GIntensity);

  Bit = 0:255;
  % FitType = 2;
  FitType = 1;

  [RFitIntensity, RExtendedBit] = FitGamma(RBit',normedRIntensity', Bit', FitType);
  % RExponent = RExtendedBit(1);
  % ROffset = RExtendedBit(2);

%   figure
%   plot(Bit,RFitIntensity)
%   hold on
%   plot(RBit,normedRIntensity, '*')
%   hold off
%   xlim([0 255])
%   ylim([0 1])
%   title('R')

  [GFitIntensity, GExtendedBit] = FitGamma(GBit',normedGIntensity', Bit', FitType);
  % GExponent = GExtendedBit(1);
  % GOffset = GExtendedBit(2);

%   figure
%   plot(Bit,GFitIntensity)
%   hold on
%   plot(GBit,normedGIntensity, '*')
%   hold off
%   xlim([0 255])
%   ylim([0 1])
%   title('G')

  inverted_gamma_params = struct('Channel',{}, 'Exponent',{}, 'Offset',{}, 'Parameter',{});

  inverted_gamma_params(1).Channel = 'R';
  % inverted_gamma_params(1).Exponent = RExponent;
  % inverted_gamma_params(1).Offset = ROffset;
  inverted_gamma_params(1).Parameter = RExtendedBit;

  inverted_gamma_params(2).Channel = 'G';
  % inverted_gamma_params(2).Exponent = GExponent;
  % inverted_gamma_params(2).Offset = GOffset;
  inverted_gamma_params(2).Parameter = GExtendedBit;

  save('inverted_gamma_params', 'inverted_gamma_params')

  RParameter = [GExtendedBit 0];
  GParameter = [RExtendedBit 0];

  correctedR = InvertGammaExtP(RParameter, 255, Bit / 255);
  correctedG = InvertGammaExtP(GParameter, 255, Bit / 255);

  correctedRBit = round(correctedR);
  correctedGBit = round(correctedG);

  inverted_gamma = struct('Channel',[repmat({'R'}, [256 1]); repmat({'G'}, [256 1])], 'Value',num2cell([Bit' / 255; Bit' / 255]), 'Bit',num2cell([Bit'; Bit']), 'CorrectedValue',num2cell([correctedRBit' / 255; correctedGBit' / 255]), 'CorrectedBit',num2cell([correctedRBit'; correctedGBit']));

  [inverted_gamma(:).Intensity] = deal([]);

  if ismember(0, RBit)
    correctedRIntensity = (RIntensity(RBit == 255) - RIntensity(RBit == 0)) / 255 * (0:255) + RIntensity(RBit == 0);
  else
    correctedRIntensity = RIntensity(RBit == 255) / 255 * (0:255);
  end

  if ismember(0, GBit)
    correctedGIntensity = (GIntensity(GBit == 255) - GIntensity(GBit == 0)) / 255 * (0:255) + GIntensity(GBit == 0);
  else
    correctedGIntensity = GIntensity(GBit == 255) / 255 * (0:255);
  end
  
  correctedIntensity = num2cell([correctedRIntensity'; correctedGIntensity']);
  [inverted_gamma.Intensity] = correctedIntensity{:};


  save('inverted_gamma', 'inverted_gamma')

end

function [Bit, normedIntensity] = NormalizeIntensity(Bit, Intensity)

  [Bit, index] = sort(Bit);
  Intensity = Intensity(index);

  if Bit(1) == 0
    normedIntensity = Intensity - Intensity(1);
  else
    normedIntensity = Intensity;
  end

  normedIntensity = normedIntensity / normedIntensity(Bit == 255);

end
