function [EIC_2D] = extractEIC(MSvaluebox,MSintbox,NbPix2ndD,MStotint,Low_value,High_value,PlotVar)

% ExtractsEIC with m/z between Low_value and High_value. PlotVar = 1
% generates a plot (if=0, then no plot displayed)

%  Low_value=337.87;
%     High_value=337.97;
    EIC=zeros(size(MSintbox,1),1);
    for k=1:length(EIC)
    	Zurb = MSintbox(k,:);
    	EIC(k) = sum(Zurb(and(MSvaluebox(k,:)>=Low_value,MSvaluebox(k,:)<=High_value)));
    end
    EIC_2D = reshape([EIC;zeros(NbPix2ndD-mod(length(MStotint),NbPix2ndD),1)],...
        NbPix2ndD,[]);
    if PlotVar
    	figure
    	pcolor(EIC_2D);
    	shading flat
    	colorbar
    end