if(~isfield(param,'ShowFig'))
    param.ShowFig=true;
end
tic
for n=1:nfnlCgIters
%     disp([Slis WhichRep n]);
    lastRes=res;
    res = fnlCg(lastRes,param);
    
    HowMuchAdvanced=grmss(lastRes-res)/ ( grmss(res) + grmss(lastRes));
%     disp(HowMuchAdvanced);
    
    im_res = param.XFM'*res;
%     figure(FigH), gmontage(abs(gflip(im_res,1))), drawnow;% title(qq)
    if(param.ShowFig)
        figure(FigH); subplot(1,2,1);gmontagef(abs(gflip(im_res,1)),'Size',[1 size(im_res,3)]); drawnow;% title(qq)
        cx=caxis;
        caxis(cx/RunFnlViewAmp);
        subplot(1,2,2);gmontagef(angle(gflip(im_res,1)),'Size',[1 size(im_res,3)]); drawnow;% title(qq)
%         xlabel([Slis WhichRep n]);
    end
    
    if(isfield(param,'RelRMSAdvanceT'))
        if(HowMuchAdvanced<param.RelRMSAdvanceT)
            break;
        end
    end
end
t=toc;