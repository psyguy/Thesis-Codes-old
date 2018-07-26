%-------------------------------------------------------------------------%
%-                           Spase Petkoski                              -%
%-Theoretical Neuroscience Group--Institut de Neurosciences des Syst?mes -%
%-                        2018, AMU Marseille                            -%
%-------------------------------------------------------------------------%
% function for calculating phases and mean field at each step of the Heun
% integration for a homogeneous connectome of KM with white noise
% theta -  phases up to history t-tau_m
% th1 -  downsampled phases, kept for all t and N nodes
% chck - counter for th
% ind -  non-zero links for each node
% idxind1 - indices for delayed values +1, for phases with non-zero delays
% idx - indices for delayed values of each phase (even those with 0 delays 
% and/or weights)
% computationaly optimized function! might be difficult to comprehend.
% doesn't work if there are 0 delays!!!  for all existent links tau>0,
% otwerwise the for loops are needed calculating theta_1
%
% no division by N in the coupling, it is taken into account in K_cnctm
function [theta, th, i, thdot] = KMcnctmHfPlosCB(theta, eta, i, chck, varstrKM)
if isfield(varstrKM, 'omega')
    omega=varstrKM.omega;
else
    omega=varstrKM.mu;
end
K_cnctm=varstrKM.K_cnctm;
h=varstrKM.h;
idx=varstrKM.idx;
idx1=varstrKM.idx1;
Nds=varstrKM.Nds;
t=varstrKM.t;
cnt=varstrKM.cnt;
N=varstrKM.N;
%-------------------------------------------------------------------------
%-------  connectome ------------
%-------------------------------------------------------------------------
kn1=omega - sum(K_cnctm.*sin(repmat(theta(:,end), 1, N) - theta(idx)),2); 
thdot=kn1+eta;
%--
kn2=omega - sum(K_cnctm.*sin(repmat(theta(:,end) + h*kn1 + eta, 1, N) - theta(idx1)),2);
%
theta(:,1:end-1)=theta(:,2:end);
theta(:,    end)=theta(:,end) + h/2*(kn1+kn2) + eta;
%---
i=i+1;
%---
if chck==i/Nds,
    th=angle(mean(exp(1j*theta( :,end-Nds:end)),2));
    if ~mod(chck,cnt)
        disp(round(t(chck)));
    end
else
    th=NaN(N,1);
end;
end

