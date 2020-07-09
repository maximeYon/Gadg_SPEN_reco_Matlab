function [FinalADCMapRaw]=ADC_QB(DiffImag, bvalues)
% for i=1:size(bvalues,2)
%     if isempty(bvalues{i})
%         bvalues{i}=0;
%     end
% end
% 
% bvalues=cell2mat(bvalues);
matrix_size=size(DiffImag);
Para.np=matrix_size(1);
Para.nv=matrix_size(2);
bvalue2D=ones(2,Para.np,Para.nv);
bvalue2D(1,:,:)=bvalues(1)*ones(Para.np,Para.nv);
bvalue2D(2,:,:)=mean(bvalues(2:end))*ones(Para.np,Para.nv);

DiffImag=double(abs(DiffImag));
DiffImag=DiffImag/max(DiffImag(:));
absiField(:,:,:,1)=(DiffImag(:,:,:,1));
if size(DiffImag,4)==2
    absiField(:,:,:,2)=DiffImag(:,:,:,2);
elseif size(DiffImag,4)==5
    absiField(:,:,:,2)=(DiffImag(:,:,:,2).*DiffImag(:,:,:,3).*DiffImag(:,:,:,4).*DiffImag(:,:,:,5)).^(1/4);
else
    absiField(:,:,:,2)=(DiffImag(:,:,:,2).*DiffImag(:,:,:,3).*DiffImag(:,:,:,4)).^(1/3);
end
absiField=double(absiField);
% absiField(:,:,:,2)=mean(DiffImag(:,:,:,2:end),4);
% absiField=absiField/max(absiField(:))*100;
%estimate adc and estimatio4n error for each imaging point
% for i=1:NumDiff
%     absiField(:,:,:,i)=abs(iField(:,:,:,i))./abs(iField(:,:,:,1));
% end

% figure()
% ha = tight_subplot(1,matrix_size(3),0.001)
fh = @(x,p)  p(1)*exp(x*p(2));
errfh = @(p,x,y) sum((y(:)-fh(x(:),p)).^2);
ADC=zeros(Para.np,Para.nv);
ADCFmin=zeros(Para.np,Para.nv);
FinalADCMap=zeros(Para.np,Para.nv,matrix_size(3));
FinalADCMapRaw=zeros(Para.np,Para.nv,matrix_size(3));
for iSlice=1:matrix_size(3)
    for m=1:Para.np
        for n=1:Para.nv
            if(absiField(m,n,iSlice,1)>0.0000001)
            x=squeeze(-bvalue2D(:,m,n));
                         
            y=log(squeeze(abs(absiField(m,n,iSlice,:))));
            p=polyfit(x,y,1);
            ADC(m,n)=p(1);
%             yfit=polyval(p,x);
            %         if (abs(S(1,m,n))>0.2*max(max(squeeze(abs(S(1,:,:))))))
            %            plot(x,y,'*b',x,yfit,'-r');
            %            pause(0.5);
            %         end
%             yresid=y-yfit;
%             sumyresid=sum(yresid.^2);
%             Sumtotal = (length(yfit)-1)* var(y);
            p0 = [exp(p(2)) p(1)];
            P = fminsearch(errfh,p0,[],x,squeeze(abs(absiField(m,n,iSlice,:))));
            ADCFmin(m,n)=P(2);
%             rsq(m,n) = 1 - sumyresid/Sumtotal;
            end
        end
    end
    
    
%     FinalADCMap(:,:,iSlice)=ADC;
    FinalADCMap(:,:,iSlice)=ADCFmin;
    FinalADCMapRaw(:,:,iSlice)=ADC;
    
end
