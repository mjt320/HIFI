function s=deichmannFormula(S0,T1_s,TR_s,TI_s,TD_s,a_rad,b_rad,N,PECentre)
%formula taken from Deichmann 2000
%S0=equilibrium signal
%TR_s=spacing between readout pulses
%TI_s=delay between inversion pulse and acquisition block
%a_rad=inversion FA (not used here - formula assumes perfect inversion)
%b_rad=readout FA
%TD_s=delay between acquisitions
%N=number of readout pulses in acquisition block
%PECentre=time of centre of k-space expressed as a fraction of the readout train duration
%s=signal

tau = TR_s .* N; %duration of the acquisition block

T1Star = ( (1/T1_s) - (1./TR_s).*log(cos(b_rad)) ).^(-1);
M0Star = S0 * ((1-exp(-TR_s/T1_s)) ./ (1-exp(-TR_s./T1Star)));

A1 = M0Star .* (1 - exp(-tau./T1Star));
A2 = S0 * (1 - exp(-TD_s/T1_s));
A3 = S0 * (1 - exp(-TI_s/T1_s));
B1 = exp(-tau./T1Star);
B2 = exp(-TD_s/T1_s);
B3 = -exp(-TI_s/T1_s);

A = A3 + A2.*B3 + A1.*B2.*B3;
B = B1.*B2.*B3;

M1 = A./(1-B);

%s = abs( ( M0Star + (M1 - M0Star).*exp(-(PECentre.*tau)./T1Star) ) .*sin(b) );
s = ( M0Star + (M1 - M0Star).*exp(-(PECentre.*tau)./T1Star) ) .*sin(b_rad);

end
