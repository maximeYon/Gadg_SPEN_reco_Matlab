function out=CombineDataDim(Data1, Data2)
% Combine two dataset different size in PE regardless difference of size
% due to not equal number of EO

SD1=size(Data1,1);
SD2=size(Data2,1);
% out=nan([SD1+SD2, gsize(Data1,[2 3 4 5])]);
if SD1<SD2
    disp('Swap You two Data !!')
    return
end

DimCat=length(size(Data1))+1;


Cat=cat(DimCat,Data1(1:SD2,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:),Data2(1:SD2,:,:,:,:,:,:,:,:,:,:,:,:,:));
out=CombineDims(Cat,[1 DimCat]);
if SD1~=SD2
    out(SD1+SD2,:,:,:,:,:,:,:,:,:,:,:,:,:)=Data1(SD1,:,:,:,:,:,:,:,:,:,:,:,:,:);
end
