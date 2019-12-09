function getSensitivity_macro(parentpath)

rng; 

% Get all folder names in the 'Results' directory
d = dir(parentpath)
isub = [d(:).isdir]; %# returns logical vector
subFolders = {d(isub).name}'; 
subFolders(ismember(subFolders,{'.','..'})) = []; 

pathList = fullfile(parentpath, subFolders)

for sub = 1:length(pathList)
    
    curr_path = pathList{sub}

    files_coneresp = dir(fullfile(curr_path, 'ConeExcitationInstances', '*.mat'));
%     files_coneresp = dir(fullfile(curr_path, 'NoisyConeExcitationInstances', '*.mat'));

    for f = 1:length(files_coneresp)
        sf_list(f) = str2double(cell2mat(extractBetween(files_coneresp(f).name, 'SF_', '_contr'))); 
        contrast_list(f) = str2double(cell2mat(extractBetween(files_coneresp(f).name, '_contr_', '.mat'))); 
    end
    SFs = unique(sf_list);
    nSF = length(SFs); 
    contrasts = unique(contrast_list); 
    nContrast = length(contrasts); 
    
    result_svm = nan(nSF, nContrast); 
    for f = 1:length(files_coneresp)
        
        curr_filename = fullfile(files_coneresp(f).folder, files_coneresp(f).name); 
        load(curr_filename, 'SVMpercentCorrect'); 
        
        val_sf = cell2mat(extractBetween(files_coneresp(f).name, 'SF_', '_contr'));
        val_contrast = str2double(cell2mat(extractBetween(files_coneresp(f).name, '_contr_', '.mat'))); 
        
        result_svm(SFs==val_sf, contrasts==val_contrast) = SVMpercentCorrect; 
%         result_svm(f) = SVMpercentCorrect; % memory of SVM


%         % Computing the sesntivitiy
%         if mod(f, nContrast) == 0
%             val_exp = cell2mat(extractBetween(files_coneresp(f).name, 'coneExcitation_noiseOff_', '_SF'));
% %             val_exp = cell2mat(extractBetween(files_coneresp(f).name, 'coneExcitation_', '_SF'));
%             
%             val_sensitivity = fitPsychometricFn(contrasts, result_svm(f-(nContrast-1):f));
%             
%             fprintf('%s SF%s: %f \n', val_exp, val_sf, val_sensitivity)
%         end
        
    end
    result_svm
end



end 