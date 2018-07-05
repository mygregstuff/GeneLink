function X=gpu(X)

% if prod([size(X)])>1e0
%     'array moved to HD due to size'
%     X=tall(X);
% else
%     
    if  gpuDeviceCount>0
    'array moved to gpu'
    X=gpuArray(X);
else
    %'no gpu detected'
end