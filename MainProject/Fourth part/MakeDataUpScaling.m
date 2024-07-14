clear; close all; clc;
%% Make data Upscaling

xxt = [];
yyt = [];
ppt = [];

for i = [32, 64, 128] 
    load(['CNN' num2str(i) '.mat'])
%     UpScale
    xt = GridResize(xt); 
    xxt = cat(4,xxt,xt) ; % concatenates(128,130,1,3*5600) 4D , All 
    yyt = [yyt ; yt] ;
    ppt = [ppt ; pt] ; 
    
end
clearvars -expect xxt ppt yyt
save('CNN_upscale.mat', '-v7.3')


function output = GridResize(xt)
    % reference number of grid for upscaling
    RNG = 128;
    for indk = 1:size(xt,4)% size(xt, 4D)
        S = xt(:,:,:,indk);
        % slicing 1st column, end column
        S = S(:, 2:end-1);
        c = RNG/length(S);
        B = nan(RNG);
        % Reshape
        for ind = 1:size(S,1) % size(s, 1(row or D-row))
            for jind = 1:size(S,2) % size(s, 2(column or D-column))
                B(c*ind-(c-1):c*ind, c*jind-(c-1):c*jind) = repmat(S(ind,jind),c,c); % y = 4x-3
            end
        end
        k = [ones(RNG,1) B ones(RNG,1)];
        output(:,:,1,indk) = k ; % make 4D 
    end

end