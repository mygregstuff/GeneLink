function X=gpu(X)
if gpuDeviceCount>0
    X=gpuArray(X);
else
    %'no gpu detected'
    X=X;
end