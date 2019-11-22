function [contrastSensitivity, contrastThreshold] = fitPsychometricFn(contrasts, dat)

%% *Introductory script illustrating how to fit a Weibul function to psychometric data (performance as a function of stimulus contrast)*
% _This tutorial fits a cumulative Weibull function to a psychometric dataset 
% using the Palamedes toolbox._
% 
% _Copyright: Nicolas P. Cottaris,  ISETBio Team, 2019_
%% 
% 
%% *Step 1.* Load an example dataset
% Load some exemplar psychometric data

% Load exemplar data
% [contrasts, fractionCorrect, nTrials] = loadPerformanceData('ExemplarPsychometricFunction');

fractionCorrect = dat./100; 
nTrials = 32; 

% % Plot the data
% plotPsychometricFunction('raw data', contrasts, fractionCorrect, [],[], [], []);
%% *Step 2a.* Set up params and options for Palamedes fitting

% Set up psychometric function model. Here we use a cumulative Weibull function
psychometricFunctionModel = @PAL_Weibull;

% Set up search grid
gridLevels = 100;
searchGridParams.alpha = logspace(log10(min(contrasts)),log10(max(contrasts)),gridLevels);
searchGridParams.beta = 10.^linspace(-4,4,gridLevels);
searchGridParams.gamma = 0.5;
searchGridParams.lambda = 0.0;

% Optimization settings for the fit
optionsParams             = optimset('fminsearch');
optionsParams.TolFun      = 1e-09;
optionsParams.MaxFunEvals = 1000;
optionsParams.MaxIter     = 1000;
optionsParams.Display     = 'off';

% Parameters for the curve fitting
% Parameters that are allowed to vary
% The parameters are: threshold, slope, guess-rate, lapse-rate
paramsFree = [1 1 0 0];
trialsNumCorrectPerContrastLevel = round(nTrials*fractionCorrect);
trialsNumPerContrastLevel = repmat(nTrials,1,length(fractionCorrect));

%% *Step 2b. Fit the data to obtain the contrast threshold*

% Fit the data and get the best fit params
paramsValues = PAL_PFML_Fit(contrasts(:), trialsNumCorrectPerContrastLevel(:), trialsNumPerContrastLevel(:), ...
            searchGridParams, paramsFree, psychometricFunctionModel, 'SearchOptions', optionsParams);
        
% Obtain the threshold at which performance cross a threshold performance, here 71%
performanceThreshold = 0.71;
contrastThreshold = psychometricFunctionModel(paramsValues, performanceThreshold, 'inverse');
contrastSensitivity = 1/contrastThreshold;
%% *Step 2c. Visualize fitted function*

% Obtain a high resolution version of the fitted function
hiResContrasts = searchGridParams.alpha;
hiResPerformance = PAL_Weibull(paramsValues, hiResContrasts);

% Plot data and fitted function
plotPsychometricFunction('fitted data', contrasts, fractionCorrect, hiResContrasts, hiResPerformance, contrastThreshold, performanceThreshold);
%% 

end 

% function [contrasts, fractionCorrect, nTrials] = loadPerformanceData(filename)
%     rootDir = tbLocateProject('ISETBioLiveScript');
%     dataFile = fullfile(rootDir, 'toolbox', 'resources', filename);
%     load(dataFile, 'contrasts', 'fractionCorrect');
%     % Assume those data were collected using 1000 trials
%     nTrials = 1000;
% end

function plotPsychometricFunction(theTitle, contrasts, fractionCorrect, hiResContrasts, hiResPerformance, contrastThreshold, performanceThreshold)
    figure(); clf;
    
    if (~isempty(hiResContrasts))
        plot(hiResContrasts, hiResPerformance, 'r-', 'LineWidth', 1.5);
        hold on
    end
    if (~isempty(contrastThreshold))
        plot(contrastThreshold*[1 1], [0 performanceThreshold], 'b-', 'LineWidth', 1.5);
        plot(contrastThreshold*[0.01 1], performanceThreshold*[1 1], 'b-', 'LineWidth', 1.5);
    end
    
    plot(contrasts,fractionCorrect, 'ko', 'MarkerSize', 12, ...
        'MarkerFaceColor', [0.7 0.5 0.5], 'MarkerEdgeColor', [0.5 0 0]);
    
    contrastTicks = [0.001 0.003 0.01 0.03 0.1 0.3 1];
    contrastTickLabels = {'.001', '.003', '.01', '.03', '.1', '.3', '1.'};
    set(gca, 'XTick', contrastTicks, 'XTickLabel', contrastTickLabels, 'YTick', 0:0.1:1.0);
    set(gca, 'XLim', [0.001 1], 'YLim', [0.4 1.05], 'XScale', 'log')
    set(gca, 'FontSize', 16)
    xlabel('\it contrast');
    ylabel('\it fraction correct');
    title(theTitle)
    axis 'square'
    grid on; box on 
end