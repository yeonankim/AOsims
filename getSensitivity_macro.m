function getSensitivity_macro(parentpath, foldername)

rng; 

% Get all folder names in the 'Results' directory
d = dir(parentpath);
isub = [d(:).isdir]; %# returns logical vector
subFolders = {d(isub).name}';
subFolders(ismember(subFolders,{'.','..'})) = []; 

pathList = fullfile(parentpath, subFolders);

fig_cnt = 0; 
for sub = 1:length(pathList)
    
    fprintf('MOSAIC CONDITION %d \n', sub);
    curr_path = pathList{sub};
    files_coneresp = dir(fullfile(curr_path, foldername, '*.mat'));

    for f = 1:length(files_coneresp)
        oi_list (f) = str2double(cell2mat(extractBetween(files_coneresp(f).name, '_oi', '_exp'))); 
        exp_list(f) = str2double(cell2mat(extractBetween(files_coneresp(f).name, '_exp', '_SF'))); 
        sf_list(f) = str2double(cell2mat(extractBetween(files_coneresp(f).name, 'SF_', '_contr'))); 
        contrast_list(f) = str2double(cell2mat(extractBetween(files_coneresp(f).name, '_contr_', '.mat'))); 
    end
    OIs = unique(oi_list);             nOI = length(OIs);
    EXPs = unique(exp_list);           nEXP = length(EXPs);
    SFs = unique(sf_list);             nSF = length(SFs); 
    contrasts = unique(contrast_list); nContrast = length(contrasts); 
    
    result_svm = nan(nSF, nContrast, nEXP, nOI); 
    for f = 1:length(files_coneresp)
        curr_filename = fullfile(files_coneresp(f).folder, files_coneresp(f).name);
        load(curr_filename, 'SVMpercentCorrect'); 
        
        val_oi = str2double(cell2mat(extractBetween(files_coneresp(f).name, '_oi', '_exp')));
        val_exp = str2double(cell2mat(extractBetween(files_coneresp(f).name, '_exp', '_SF')));
        val_sf = str2double(cell2mat(extractBetween(files_coneresp(f).name, 'SF_', '_contr'))); 
        val_contrast = str2double(cell2mat(extractBetween(files_coneresp(f).name, '_contr_', '.mat'))); 
        
        result_svm(SFs==val_sf, contrasts==val_contrast, EXPs==val_exp, OIs==val_oi) = SVMpercentCorrect; 
    end
    
    for k = 1:nOI
        for i = 1:nEXP
            for j = 1:nSF
                fig_cnt = fig_cnt + 1; 
                col_valid = ~isnan(result_svm(j,:,i, k));
                [val_sensitivity, cont_thres] = fitPsychometricFn(contrasts(col_valid), result_svm(j,col_valid,i,k));
                fprintf('%d. OI%d EXP%d SF%f: %f (%f) \n', fig_cnt, OIs(k), EXPs(i), SFs(j), val_sensitivity, cont_thres);
            end
        end
    end

end

end 