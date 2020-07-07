function Out=SmoothBySlices(In,Params,Sigma)
% function Out=SmoothBySlices(In,Params,Sigma)

h=fspecial('gaussian',Params,Sigma);
Out=In;
for i=1:size(In,3)
    for j=1:size(In,4)
        for k=1:size(In,5)
            for l=1:size(In,6)
                for m=1:size(In,7)
                    for n=1:size(In,8)
                        CurI=squeeze(In(:,:,i,j,k,l,m,n));
                        Out(:,:,i,j,k,l,m,n)=imfilter(CurI, h);
                    end
                end
            end
        end
    end
end