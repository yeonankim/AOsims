function getSensitivity_macro(parentpath, foldername, displayfilename, varargin)
%varargin to change the order of mosaic to plot, or to analyze some of the
%mosiac patterns. 

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
mosaic_compute_num = all_mosaic_num;
mosaic_compute_order = 1:length(pathList);

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
        mosaic_cond_name{sub} = sprintf('L:M = %g:%g', round(Ldensity/min(Ldensity, Mdensity),2), round(Mdensity/min(Ldensity, Mdensity),2)); 
        fprintf('MOSAIC CONDITION %d (%s) \n', mosaic_compute_num(sub), mosaic_cond_name{sub});
        
        
        % Make the list of files for the given simulation condition
        files_coneresp = dir(fullfile(curr_path, foldername, '*.mat'));
        
        % Go through file names and extract the # layers of OI, EXP, SF, and CONTRAST conditions
        condmatrix = nan(length(files_coneresp), 6); 
        for f = 1:length(files_coneresp)
            val_oi = str2double(cell2mat(extractBetween(files_coneresp(f).name, '_oi', '_exp')));
            val_exp = str2double(cell2mat(extractBetween(files_coneresp(f).name, '_exp', '_SF')));
            val_sf = str2double(cell2mat(extractBetween(files_coneresp(f).name, 'SF_', '_contr')));
            val_contrast = str2double(cell2mat(extractBetween(files_coneresp(f).name, '_contr_', '.mat')));
            
            curr_filename = fullfile(files_coneresp(f).folder, files_coneresp(f).name);
            load(curr_filename, 'SVMpercentCorrect');
            nSVMrep = length(SVMpercentCorrect);
            
            condmatrix(f,:) = [val_oi, val_exp, val_sf, val_contrast, nSVMrep, f]; 
            svmresult{f} = SVMpercentCorrect;
        end
        OIs = unique(condmatrix(:,1));             nOI = length(OIs);
        EXPs = unique(condmatrix(:,2));           nEXP = length(EXPs);
        SFs = unique(condmatrix(:,3));             nSF = length(SFs);
        contrasts = unique(condmatrix(:,4)); nContrast = length(contrasts);
        nSVM = min(condmatrix(:,5)); % NEED TO CHECK MAX?
        
        for k = 1:nOI
            for i = 1:nEXP
                for j = 1:nSF
                    fig_cnt = fig_cnt + 1; 
                    this_ = condmatrix(condmatrix(:,1)==OIs(k) & condmatrix(:,2)==EXPs(i) & condmatrix(:,3)==SFs(j), [4,6]);
                    this_contrast = contrasts(ismember(contrasts, this_(:,1)));
                    this_svm = reshape(cell2mat(svmresult(this_(:,2))), [nSVM, length(this_contrast)]);
                    
                    if nSVM > 1
                        sensitivity = zeros(nSVM,1);
                        thres = zeros(nSVM,1);
                        for l = 1:nSVM
                            [sensitivity(l), thres(l)] = fitPsychometricFn(this_contrast, this_svm(l,:));
                        end
                        val_sensitivity = mean(sensitivity);
%                         val_sensitivity_se = std(sensitivity)./sqrt(nSVM);
                        cont_thres = mean(thres);
                        
%                         res_sensitivity(j,sub,i,k) = mean(sensitivity); 
                        res_sensitivity_se(j,sub,i,k) = std(sensitivity)./sqrt(nSVM);
                    else
                        [val_sensitivity, cont_thres] = fitPsychometricFn(this_contrast, this_svm);
%                         res_sensitivity(j,sub,i,k) = val_sensitivity; 
                    end
                    [L_cont, M_cont] = computeConeContrast(displayfilename, SFs(j), cont_thres);
                    fprintf('%d. OI%d EXP%d SF%d: %f (michelson %f, L cone %f, M cone %f) \n', fig_cnt, OIs(k), EXPs(i), SFs(j), val_sensitivity, cont_thres, L_cont, M_cont);
                    res_sensitivity(j,sub,i,k) = val_sensitivity; 
                    res_Lcont(j,sub,i,k) = L_cont;
                    res_Mcont(j,sub,i,k) = M_cont;
                end
            end
        end
        
    else
        fprintf('%d does not exist here: Skipping. \n', foldername);
        
    end %if result folder exists

end

% save('result_save_temp.mat', 'mosaic_cond_name', 'mosaic_compute_num', 'SFs', 'res_sensitivity', 'res_Lcont', 'res_Mcont');
if exist('res_sensitivity_se', 'var')
%     save('result_save_temp.mat', 'res_sensitivity_se', '-append');
    plotSensitivity(res_sensitivity, mosaic_cond_name, mosaic_compute_num, SFs, res_sensitivity_se); 
else
    plotSensitivity(res_sensitivity, mosaic_cond_name, mosaic_compute_num, SFs); 
end

end 


function plotSensitivity(result, mosaic_cond_name, mosaic_compute_num, SFs, varargin)
%% Plot the results 
% such that each mosaic condition gets assigned a consistent color
plotcolors = lines; 
mosaic_cond_name(cellfun('isempty',mosaic_cond_name)) = [];
for k = 1:size(result,4)
    figure('Name',sprintf('Optics condition %d', k),'NumberTitle','off'); %sensitivity
    for i = 2:size(result,3)
%        subplot(1,3,i); 
%        b=bar(result(:,:,i,k)); 
%        for m = 1:length(mosaic_compute_num); b(m).FaceColor = plotcolors(mosaic_compute_num(m),:); end

        subplot(1,2,i-1); 
        if nargin == 4
            p = plot(result(:,:,i,k), '.-', 'LineWidth', 1, 'MarkerSize', 12);
        elseif nargin == 5
            p = errorbar(result(:,:,i,k), varargin{1}(:,:,i,k), '.-', 'LineWidth', 1, 'MarkerSize', 12);
        end
        for m= 1:length(mosaic_compute_num); p(m).Color = plotcolors(mosaic_compute_num(m),:); end; xticks(1:length(SFs)); 
        
        xticklabels(string(SFs)); title(sprintf('EXP%d', i-1)); if i == 3; L = legend(mosaic_cond_name); L.Title.String = 'L:M ratio'; end
    end
end 
end