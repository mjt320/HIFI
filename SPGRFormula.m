function s=SPGRFormula(S0,T1_s,TR_s,b_rad)
%standard formula for steady-state spoiled gradient echo signal
%S0=equilibrium signal
%T1=T1 relaxation time
%TR=spacing between readout pulses
%b=readout FA

s=S0 .* (((1-exp(-TR_s./T1_s)).*sin(b_rad)) ./ (1-exp(-TR_s./T1_s).*cos(b_rad)) );

end
