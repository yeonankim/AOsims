function [L_cont, M_cont] = computeConeContrast(displayfilename, sf, cont)

presentationDisplay = displayCreate(displayfilename);
[scene, I] = generateGaborSceneAO(presentationDisplay, 1, 1, sf, cont); 
% figure; hold on; 
% plot(I(round(size(I,1)/2),:, 1), 'r');
% plot(I(round(size(I,1)/2),:, 2), 'g');
% 
% figure; hold on; 
% plot(linearizedI(round(size(linearizedI,1)/2),:, 1), 'r');
% plot(linearizedI(round(size(linearizedI,1)/2),:, 2), 'g');

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
if max(unique(map_l(loc_r))) >= max(unique(map_l(loc_g)))
    L_r = max(unique(map_l(loc_r))); 
    L_g = min(unique(map_l(loc_g))); 
    L_cont = (L_r - L_g)/(L_r + L_g);
else
    L_r = min(unique(map_l(loc_r)));
    L_g = max(unique(map_l(loc_g))); 
    L_cont = (L_g - L_r)/(L_r + L_g);
end
if max(unique(map_m(loc_g))) >= max(unique(map_m(loc_r)))
    M_g = max(unique(map_m(loc_g))); 
    M_r = min(unique(map_m(loc_r))); 
    M_cont = (M_g - M_r)/(M_g + M_r);
else
    M_g = min(unique(map_m(loc_g)));
    M_r = max(unique(map_m(loc_r))); 
    M_cont = (M_r - M_g)/(M_g + M_r);
end

end