function [L_cont, M_cont, L_rms, M_rms] = computeConeContrast(displayfilename, sf, cont)

presentationDisplay = displayCreate(displayfilename);
[scene, I] = generateGaborSceneAO(presentationDisplay, 1, 1, sf, cont); 
% [scene, I] = generateGaborSceneAO(presentationDisplay, 0, 1, sf, cont); 

% figure; hold on; 
% plot(I(round(size(I,1)/2),:, 1), 'r');
% plot(I(round(size(I,1)/2),:, 2), 'g');

map_lms = sceneGet(scene, 'lms');

% figure; hold on; 
% plot(map_lms(round(size(map_lms,1)/2),:, 1), 'r');
% plot(map_lms(round(size(map_lms,1)/2),:, 2), 'g');

%% ---- do above meanline instead
loc_r = I(:,:,1) == max(max(I(:,:,1)));
loc_g = I(:,:,2) == max(max(I(:,:,2))); 

%% ---- do above meanline instead
map_l = map_lms(:,:,1);
map_m = map_lms(:,:,2);

Lmax = max(max(map_lms(:,:,1)));
Lmin = min(min(map_lms(:,:,1)));
L_cont = (Lmax-Lmin)./(Lmax+Lmin);

Mmax = max(max(map_lms(:,:,2)));
Mmin = min(min(map_lms(:,:,2)));
M_cont = (Mmax-Mmin)./(Mmax+Mmin);

% if max(unique(map_l(loc_r))) >= max(unique(map_l(loc_g)))
%     L_r = max(unique(map_l(loc_r))); 
%     L_g = min(unique(map_l(loc_g))); 
%     L_cont = (L_r - L_g)/(L_r + L_g);
% else
%     L_r = min(unique(map_l(loc_r)));
%     L_g = max(unique(map_l(loc_g))); 
%     L_cont = (L_g - L_r)/(L_r + L_g);
% end
% if max(unique(map_m(loc_g))) >= max(unique(map_m(loc_r)))
%     M_g = max(unique(map_m(loc_g))); 
%     M_r = min(unique(map_m(loc_r))); 
%     M_cont = (M_g - M_r)/(M_g + M_r);
% else
%     M_g = min(unique(map_m(loc_g)));
%     M_r = max(unique(map_m(loc_r))); 
%     M_cont = (M_r - M_g)/(M_g + M_r);
% end
L_rms = sqrt((1/numel(map_l))*sum(sum((map_l-mean(map_l(:))).^2)));
M_rms = sqrt((1/numel(map_m))*sum(sum((map_m-mean(map_m(:))).^2)));

end