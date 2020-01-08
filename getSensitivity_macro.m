function getSensitivity_macro(parentpath, foldername, displayfilename, varargin)

rng; 

% Get all mosaicCond folder names in the 'Results' directory
d = dir(parentpath);
isub = [d(:).isdir]; %# returns logical vector
subFolders = {d(isub).name}';
subFolders(ismember(subFolders,{'.','..'})) = []; 

pathList = fullfile(parentpath, subFolders);
for sub = 1:length(pathList)
    all_mosaic_num(sub) = str2double(extractAfter(pathList{sub}, 'mosaicCond'));
end
if nargin == 4
    mosaic_compute_num = varargin{1};  % Contains mosaicCond number to compute
    mosaic_compute_order = [];         % Index to pahtList for each mosaicCond number in 'mosaic_custom_order'
    for o = 1:length(mosaic_compute_num)
        id = find(all_mosaic_num == mosaic_compute_num(o));
        mosaic_compute_order = [mosaic_compute_order, id];
        if isempty(id)
            mosaic_compute_num = mosaic_compute_num(mosaic_compute_num~=id); 
            fprintf('Warning: mosaicCond%d does not exist, Skipping.', mosaic_custom_order(o));
        end
    end
elseif nargin < 4
    mosaic_compute_num = all_mosaic_num;
    mosaic_compute_order = 1:length(pathList);
else 
    error('Too many input argument'); 
end 

fig_cnt = 0;
for sub = 1:length(mosaic_compute_order) % For each mosaic
    
    % For a give mosaic foler
    curr_path = pathList{mosaic_compute_order(sub)};
    
    if exist(fullfile(curr_path, foldername), 'dir') == 7 % Only if the given simulation condition exsit for this mosaic
        
        % Extract the L:M ratio for this mosaic and print on command line
        mosaicfile = dir(fullfile(curr_path, '*.mat')); 
        Ldensity = str2double(cell2mat(extractBetween(mosaicfile.name, 'L', 'M')));
        Mdensity = str2double(cell2mat(extractBetween(mosaicfile.name, 'M', 'S')));
        fprintf('MOSAIC CONDITION %d (L:M = %d:%d) \n', mosaic_compute_num(sub), Ldensity/min(Ldensity, Mdensity), Mdensity/min(Ldensity, Mdensity));
        
        % Make the list of files for the given simulation condition
        files_coneresp = dir(fullfile(curr_path, foldername, '*.mat'));
        
        % Go through file names and extract the # layers of OI, EXP, SF, and CONTRAST conditions
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
        
        % The matrix to which SVM result will be assigned 
        result_svm = nan(nSF, nContrast, nEXP, nOI);
        % Go through each file again, this time to extract SVM results
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
                    [L_cont, M_cont] = computeConeContrast(displayfilename, SFs(j), cont_thres);
                    fprintf('%d. OI%d EXP%d SF%d: %f (michelson %f, L cone %f, M cone %f) \n', fig_cnt, OIs(k), EXPs(i), SFs(j), val_sensitivity, cont_thres, L_cont, M_cont);
                    
                    res_sensitivity(j, sub, i, k) = val_sensitivity; 
                    res_Lcont(j, sub, i, k) = L_cont;
                    res_Mcont(j, sub, i, k) = M_cont;
                end
            end
        end
        
    end %if result folder exists

end


%% Plot the results 
% such that each mosaic condition gets assigned a consistent color
barcolors = lines; 
for k = 1:nOI
    f1 = figure; %sensitivity
    f2 = figure; %L-cone contrast
    f3 = figure; %M-cone contrast
    for i = 1:nEXP
        figure(f1); subplot(1,3,i); b1=bar(res_sensitivity(:,:,i,k)); xticklabels({'SF10', 'SF30', 'SF60'}); title(sprintf('EXP%d', i));
        figure(f2); subplot(1,3,i); b2=bar(res_Lcont(:,:,i,k)); xticklabels({'SF10', 'SF30', 'SF60'}); title(sprintf('EXP%d', i));
        figure(f3); subplot(1,3,i); b3=bar(res_Mcont(:,:,i,k)); xticklabels({'SF10', 'SF30', 'SF60'}); title(sprintf('EXP%d', i));
        for m = 1:length(mosaic_compute_num)
            b1(m).FaceColor = barcolors(mosaic_compute_num(m),:);
            b2(m).FaceColor = barcolors(mosaic_compute_num(m),:);
            b3(m).FaceColor = barcolors(mosaic_compute_num(m),:);
        end
    end
end 