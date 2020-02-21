function SensMaps = RunESPIRiTForSensMaps(SensFull,ncalib,ksize)

DATA=ifft2c(SensFull);
[sx,sy,Nc] = size(DATA);
% ncalib = 40; % use 24 calibration lines to compute compression
% ksize = [6,6]; % kernel size


% Threshold for picking singular vercors of the calibration matrix
% (relative to largest singlular value.

eigThresh_1 = 0.02;

% threshold of eigen vector decomposition in image space.
eigThresh_2 = 0.9;

% crop a calibration area
calib = crop(DATA,[ncalib,ncalib,Nc]);


[k,S] = dat2Kernel(calib,ksize);
idx = max(find(S >= S(1)*eigThresh_1));

[M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);

SensMaps = M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_2,[1,1,Nc]);
SensMaps=rot90(SensMaps,2);

end