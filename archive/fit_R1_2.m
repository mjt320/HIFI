function [T1,S0,k,modelFit,R1_LCI,R1_UCI,RSq,exitFlag]=fit_R1_2(S,isFitIR,isFit,TR_s,FA_rad,TI_s,PECentre,NReadout,NTry,options)
%fit T1

tic;

if exist('options','var')
    tolFun=options.tolFun;
    tolX=options.tolX;
else
    tolFun=1e-6;
    tolX=1e-6;
end

y=S(isFit); %select only parts of signal to be fitted

T1=nan; S0=nan; modelFit=nan(1,sum(isFit)); k=nan; RSq=nan; model=nan(1,sum(isFit)); R1_LCI=nan; R1_UCI=nan; exitFlag=nan; %initialise variables

isFitSPGR=~logical(isFitIR) & logical(isFit); %images that are SPGR and should be fitted
if sum(isFitIR)>0; isHIFI=1; else isHIFI=0; end %if there are no IR- scans, this is variable flip angle (VFA) not HIFI (and vice versa)
NScans=size(isFit,2);
NFitIRScans=sum(isFitIR);

x0=[nan nan]; xLower=[0 0]; xUpper=[inf inf]; %set initial parameters and constraints (T1, S0, k)
if isHIFI; x0(3)=1; xLower(3)=0; xUpper(3)=inf; end

%%nested function to generate signal based on parameters, called by lsqcurvefit
    function s=calcSignal(c,t) % c(1)=T1 c(2)=S0 c(3)=k
        if size(c,2)==2; c(3)=1; end
        s=nan(1,NScans);
        s(isFitSPGR)=abs(SPGRFormula(c(2),c(1),TR_s(isFitSPGR),c(3)*FA_rad(isFitSPGR))); %calculate SPGR signals
        s(isFitIR)=abs(deichmannFormula(c(2),c(1),TR_s(isFitIR),TI_s(isFitIR),zeros(NFitIRScans,1),pi*ones(NFitIRScans,1),c(3)*FA_rad(isFitIR),NReadout(isFitIR),PECentre(isFitIR))); %calculate IR-SPGR signals
        s=s(isFit); %output only the signals to be fitted
    end
signalFunction=@calcSignal;

%quickly estimate initial parameters based on first and last SPGR scans
n1=min(find(isFitSPGR)); n2=max(find(isFitSPGR)); TR_2FA=TR_s(n1); a1=FA_rad(n1); a2=FA_rad(n2);
S1=S(n1); S2=S(n2);
SR=S1/S2;
x0(1)=TR_2FA/log((SR*sin(a2)*cos(a1) - sin(a1)*cos(a2))/(SR*sin(a2) - sin(a1))); %initial T1 guess
x0(2)=S1 * ( (1-exp(-TR_2FA/x0(1))*cos(a1)) / ((1-exp(-TR_2FA/x0(1)))*sin(a1)) ); %initial S0 guess

%check for invalid starting values
if x0(1)<=0 || imag(x0(1))~=0 || isnan(x0(1)) || isinf(x0(1))...
        || x0(2)<=0 || imag(x0(2))~=0 || isnan(x0(2)) || isinf(x0(2));
    x0(1)=1; x0(2)=1000;
end

%% do non-linear fitting
RSqTry=nan(1,NTry); out=cell(1,NTry); exitFlags=nan(1,NTry); jac=cell(1,NTry); resid=cell(1,NTry);
if isHIFI x=nan(NTry,3); else x=nan(NTry,2); end
for iTry=1:NTry %loop through fitting attempts
    if iTry>1; x0_final=x0.*(1.5*rand(size(x0))+0.5); else x0_final=x0; end; %use previous estimates as initial guess for first attempt, then vary initial parameters randomly
    [x(iTry,:),resnorm,resid{iTry},exitFlags(iTry),out{iTry},lambda,jac{iTry}]=lsqcurvefit(signalFunction,x0_final,[],y... %call lsqcurvefit
        ,xLower,xUpper,optimset('Display','off','TypicalX',x0,'TolFun',tolFun,'TolX',tolX));
    
    RSqTry(iTry)=1 - sum((y-squeeze(calcSignal(x(iTry,:)))).^2) / (sum((y-mean(y)).^2)); %calculate RSq
end
[RSq,bestIdx]=max(RSqTry); %find best fitting attempt based on RSq

if exitFlags(bestIdx)>0 %assign best results to output variables
    T1=x(bestIdx,1); S0=x(bestIdx,2);
    if isHIFI; k=x(bestIdx,3); end
    modelFit=calcSignal(x(bestIdx,:));
    ci = nlparci(x(bestIdx,:),resid{bestIdx},'jacobian',jac{bestIdx}); %confidence intervals
    R1_LCI=1/ci(1,2); R1_UCI=1/ci(1,1);
    exitFlag=exitFlags(bestIdx);
end

timeElapsed=toc;

if rand<0.01 %randomly plot data to check it's working
    figure(1),plot(1:sum(isFit),y,'ko',1:sum(isFit),modelFit,'b-'); xlim([0 NScans+1]);
    title({['Initial coefficients: ' num2str(x0)] ['Fitted coefficients: ' num2str(x(bestIdx,:))] ['Time elapsed: ' num2str(timeElapsed)] ['Exit flag: ' num2str(exitFlag) ' RSq: ' num2str(RSq)] ['Evaluations: ' num2str(out{bestIdx}.funcCount)]});
    pause(0.1);
end

end
