function Out=CombineDims(In,Dims)
if(size(In,Dims(2))==1)
    Out=In;
    return;
end


MinDim=min(Dims);
AllDims=1:ndims(In);
if(any(gsize(In,Dims)==1))
    MaxDim=max(Dims);
    Out=permute(In,[1:(MinDim-1) MaxDim (MinDim+1):(MaxDim-1) MinDim (MaxDim+1):ndims(In)]);
    return;
end
if(numel(Dims)==2 && Dims(1)~=MinDim)
    Out=CombineDims(permute(In,[1:(MinDim-1) Dims(1) (MinDim+1):(Dims(1)-1) MinDim (Dims(1)+1):ndims(In)]),sort(Dims));
    return;
end
OtherDims=setdiff(AllDims,Dims);

Out=permute(In,[OtherDims Dims(end:-1:1)]);
Out=reshape(Out,[gsize(In,OtherDims) prod(gsize(In,Dims))]);
Out=permute(Out,[1:(MinDim-1) AllDims(end)-1 Dims(1):(AllDims(end)-2)]);
