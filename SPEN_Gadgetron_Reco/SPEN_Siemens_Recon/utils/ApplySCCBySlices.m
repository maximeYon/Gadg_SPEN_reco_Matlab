function Out=ApplySCCBySlices(In,Mtx,ncc)
nSlices=size(In,4);
Out=zeros([gsize(In,1:2) ncc gsize(In,4:11)]);
for s=1:nSlices
    Out(:,:,:,s,:,:)=permute(MultMatTensor(Mtx(:,1:ncc,s).',permute(In(:,:,:,s,:,:),[3 1 2 4 5 6 7 8])),[2 3 1 4 5 6 7 8]);
end