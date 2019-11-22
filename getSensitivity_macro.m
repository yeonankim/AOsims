function getSensitivity_macro(parentpath, contrasts)

rng; 
nContrast = length(contrasts);

% Get all folder names in the 'Results' directory
d = dir(parentpath);
isub = [d(:).isdir]; %# returns logical vector
subFolders = {d(isub).name}';
subFolders(ismember(subFolders,{'.','..'})) = []; 

pathList = fullfile(parentpath, subFolders); 

for sub = 1:length(pathList)
    
    curr_path = pathList{sub};
    
    files_coneresp = dir(fullfile(curr_path, 'ConeExitationInstances', '*.mat'));
    for f = 1:length(files_coneresp)
        
        curr_filename = fullfile(files_coneresp(f).folder, files_coneresp(f).name); 
        load(curr_filename);
        
        % Computing the sesntivitiy
        result_svm(f) = SVMpercentCorrect; % memory of SVM
        if mod(f, nContrast) == 0
            val_exp = extractBetween(s, 'coneExcitation_noiseOff_', '_SF');
            val_sf = extractBetween(s, 'SF_', '_contr');
            val_sensitivity = fitPsychometricFn(contrasts, svm_pcorr(f-5:f));
            
            fprintf('%s SF%s: %f \n', val_exp, val_sf, val_sensitivity)
        end
        
    end
end

end 

