function [Out, Factor]=NormalizeByRefImages(ToNormalize,Ref,Range,Sig)
RefS=SmoothBySlices(abs(Ref),Range,Sig);
ToNormalizeS=SmoothBySlices(abs(ToNormalize),Range,Sig);

Out=ToNormalize;
Factor=NaN(gsize(ToNormalize,3:4));
for i=1:size(ToNormalize,3)
    for j=1:size(ToNormalize,4)
        for k=1:size(ToNormalize,5)
            Factor(i,j,k)=Row(RefS(:,:,i,j,k))/Row(ToNormalizeS(:,:,i,j,k));
            Out(:,:,i,j,k)=ToNormalize(:,:,i,j,k)*Factor(i,j,k);
        end
    end
end