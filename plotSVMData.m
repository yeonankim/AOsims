load('mosaic_cond2.mat'); 

fontSize = 15; 

contrasts = linspace(0.001, 0.03, 6); 
sf_range = linspace(5, 50, 7); 

clear xlabel ylabel
for i = 1:size(Result, 2)
    
    
    figure; 
%     subplot(1,3,1); plot(1:6, Result{i}(1:6, :), '-o', 'LineWidth', 2); xticklabels(string(contrasts)); ylabel({'SVM performance'; '(% correct)'});  set(gca,'fontsize', fontSize);
%     subplot(1,3,2); plot(1:6, Result{i}(7:12, :), '-o', 'LineWidth', 2); xticklabels(string(contrasts)); xlabel('Spatial Frequency');           set(gca,'fontsize', fontSize);
%     subplot(1,3,3); plot(1:6, Result{i}(13:18, :), '-o', 'LineWidth', 2); xticklabels(string(contrasts)); legend(string(sf_range));             set(gca,'fontsize', fontSize);
    
    subplot(1,3,1); plot(1:6, reshape(reshape(Result{i}(1:6,:)', [1, 42]), [6,7]), '-o', 'LineWidth', 2); xticklabels(string(contrasts)); ylabel({'SVM performance'; '(% correct)'});  set(gca,'fontsize', fontSize);
    subplot(1,3,2); plot(1:6, reshape(reshape(Result{i}(7:12,:)', [1, 42]), [6,7]), '-o', 'LineWidth', 2); xticklabels(string(contrasts)); xlabel('Contrasts');           set(gca,'fontsize', fontSize);
    subplot(1,3,3); plot(1:6, reshape(reshape(Result{i}(13:18,:)', [1, 42]), [6,7]), '-o', 'LineWidth', 2); xticklabels(string(contrasts)); legend(string(sf_range));             set(gca,'fontsize', fontSize);
end 


