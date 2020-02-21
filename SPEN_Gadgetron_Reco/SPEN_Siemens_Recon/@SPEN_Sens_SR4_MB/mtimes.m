function res = mtimes(a,b)
% res = mtimes(FT, x)
%


if a.adjoint % Sig -> image
    if(a.DoROFT(1))
        resa = reshape(b,a.dataSize(1),a.DoROFT(2),a.nChannels);
%         resa =
%         zpad(resa,a.dataSize(1),a.imSize(2),size(resa,3),size(resa,4));
%         SFC MODIFICATION
        %%
        resa = zpad(resa,[a.dataSize(1),a.imSize(2),size(resa,3)]);
        %%
%         resa=padLeft(resa,a.dataSize(2)-a.DoROFT(2),2);
        resa = ifft1cg(resa,2);
    else
        resa = reshape(b,a.dataSize(1),a.dataSize(2),a.nChannels);
    end
    res=zeros(a.imSize);
    for d=1:a.nBands
        for c=1:a.nChannels
            res(:,:,c,d)=MultMatTensor(a.SRMat(:,:,c,d)',resa(:,:,c));
        end
    end
    res=squeeze(sum(permute(repmat(res,[1 1 1 1]),[1 2 4 3]).*a.invSens,4));
%     res=squeeze(sum(permute(repmat(res,[1 1 1 a.nBands]),[1 2 4 3]).*a.invSens,4));
    
    switch a.mode
        case 0
            res = real(res);
        case 1
            res = exp(1i*angle(res));
    end
  
else % image -> signal
    bb = reshape(b,a.imSize(1),a.imSize(2),a.nBands);
    
    switch a.mode
        case 0
            bb = real(bb);
        case 1
            bb = exp(1i*angle(bb));
    end
    
    bb = permute(repmat(bb,[1 1 1 a.nChannels]),[1 2 4 3]).*a.Sens; % phase correct
    for d=1:a.nBands
        for c=1:a.nChannels
            res(:,:,c,d)=MultMatTensor(a.SRMat(:,:,c,d),bb(:,:,c,d));
        end
    end
    if(a.DoROFT(1))
        res = fft1cg(res,2);
        res = cropg(res,a.dataSize(1),a.DoROFT(2));
    end
    
%     res=squeeze(sum(repmat(a.SRMat,[1 1 gsize(bb,2:3)]).*repmat(permute(bb,[4 1 2 3]),[size(bb,1) 1 1 1]),2));
%     res=squeeze(sum(repmat(a.SRMat,[1 1 gsize(bb,2:3)]).*repmat(permute(bb,[4 1 2 3]),[size(a.SRMat,1) 1 1 1]),2));
    
    
%     res=MultMatTensor(diag(a.SRMat(100,:)),bb);
%     res=fft1cg(res,1);
%     res=res(1:a.PEFactor:end,:,:,:);
    res=sum(res,4);
end

if size(b,2) == 1
    res = res(:);
end



    
