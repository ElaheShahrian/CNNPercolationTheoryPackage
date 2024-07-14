clear; close all; clc;
%% Merge all
xt = [] ;
yt = [] ;
pt = [] ;
for k = [128] 
    for i = [30, 40 , 50 , 60, 70 , 80, 90]
        load(['p' num2str(i) '_Data_' num2str(k) '.mat'])
        xt = cat(4,xt,XT) ; % concatenates(128, 130, 1, 800(+800(*7))) 4D
        yt = [yt ; YT] ;
        pt = [pt ; PT] ;
    end
end

clearvars -except pt yt xt
save("CNN128.mat")
