function getSensitivity_macro(parentpath)

rng; 

% Get all folder names in the 'Results' directory
d = dir(parentpath);
isub = [d(:).isdir]; %# returns logical vector
subFolders = {d(isub).name}'; 
subFolders(ismember(subFolders,{'.','..'})) = []; 

pathList = fullfile(parentpath, subFolders);

for sub = 1:length(pathList)
    
    curr_path = pathList{sub};

    files_coneresp = dir(fullfile(curr_path, 'ConeExcitationInstances', '*.mat'));
%     files_coneresp = dir(fullfile(curr_path, 'NoisyConeExcitationInstances', '*.mat'));

    for f = 1:length(files_coneresp)
        exp_list(f) = str2double(cell2mat(extractBetween(files_coneresp(f).name, 'coneExcitation_noiseOff_exp', '_SF'))); 
        sf_list(f) = str2double(cell2mat(extractBetween(files_coneresp(f).name, 'SF_', '_contr'))); 
        contrast_list(f) = str2double(cell2mat(extractBetween(files_coneresp(f).name, '_contr_', '.mat'))); 
    end
    EXPs = unique(exp_list);           nEXP = length(EXPs);
    SFs = unique(sf_list);             nSF = length(SFs); 
    contrasts = unique(contrast_list); nContrast = length(contrasts); 
    
    result_svm = nan(nSF, nContrast, nEXP); 
    for f = 1:length(files_coneresp)
        
        curr_filename = fullfile(files_coneresp(f).folder, files_coneresp(f).name); 
        load(curr_filename, 'SVMpercentCorrect'); 
        
        val_exp = str2double(cell2mat(extractBetween(files_coneresp(f).name, 'coneExcitation_noiseOff_exp', '_SF')));
%         val_exp = str2double(cell2mat(extractBetween(files_coneresp(f).name, 'coneExcitation_exp', '_SF')));
        val_sf = str2double(cell2mat(extractBetween(files_coneresp(f).name, 'SF_', '_contr'))); 
        val_contrast = str2double(cell2mat(extractBetween(files_coneresp(f).name, '_contr_', '.mat'))); 
        
        result_svm(SFs==val_sf, contrasts==val_contrast, EXPs==val_exp) = SVMpercentCorrect; 
    end
    
    for i = 1:nEXP
        for j = 1:nSF
            val_sensitivity = fitPsychometricFn(contrasts(~isnan(result_svm(j,:,i))), result_svm(j,:,i));
            fprintf('EXP%d SF%f: %f \n', EXPs(i), SFs(j), val_sensitivity); 
        end
    end

end

end 