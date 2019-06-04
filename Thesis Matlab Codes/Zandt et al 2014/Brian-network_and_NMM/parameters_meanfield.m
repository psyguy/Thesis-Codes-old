function [p] = parameters_meanfield(p,W2,W3,W4,W5)
%Generate a structure p that contains all parameters, that is called with p.parameter

%% simulation parameters
p.Dt = 1e-3; % consider 0.1e-3 for better accuracy; only used for stepping algorithm
p.T = [0,2];  %[s] time axis of simulation
p.maxstep = [];
p.solvewithode23s = false;  % stepping or ode23s algorithm? (ode23 cannot handle noise)

%% synaptic conductances
Nnm = 15;
p.h0 = zeros(1,Nnm);
%p.h0(1) = 0;      % [mS/cm^2], initial synaptic conductance
%p.h0(2) = 0;      % [mS/cm^2], initial synaptic conductance derivative
%p.h0(3) = 0;      % [mS/cm^2], initial synaptic conductance std

p.Vth = -55; p.Eex = 50; p.Ein = -82;
p.Cinh = (p.Ein - p.Vth)/(p.Eex - p.Vth); %effective conductance of inhibitory synapses in terms of excitatory synapse
%p.Cinh = -0.238; 

%% Calculate constants for approximating H with an exponential
% used to calculate sigma_g
p.g0_approx(1) = exp(1)/sqrt(6)*p.g0(1);
tau = 1./p.gamma;
p.tau_approx(1) =3*tau(1);

p.g0_approx(2:5) = exp(1)/2*p.g0(2:5);
p.tau_approx(2:5) = 2*tau(2:5);


%% connectivity parameters
% determine macroscopip connectivity parameters from the connectivity
% matrices generated in Brian

p.Nsyn(1) = NaN;
p.Nsyn(2) = mean(sum(W2,2));    % E -> E
p.Nsyn(3) = mean(sum(W3,2));    % E -> I
p.Nsyn(4) = mean(sum(W4,2));    % I -> E
p.Nsyn(5) = mean(sum(W5,2));    % I -> I

p.Nsyna(1) = NaN;
p.Nsyna(2) = sum(var(W2-repmat(mean(W2,2),1,size(W2,2)),1,1));
p.Nsyna(3) = sum(var(W3-repmat(mean(W3,2),1,size(W3,2)),1,1));
p.Nsyna(4) = sum(var(W4-repmat(mean(W4,2),1,size(W4,2)),1,1));
p.Nsyna(5) = sum(var(W5-repmat(mean(W5,2),1,size(W5,2)),1,1));

p.varN(1) = NaN;
p.varN(2) = var(sum(W2,2),1);
p.varN(3) = var(sum(W3,2),1);
p.varN(4) = var(sum(W4,2),1);
p.varN(5) = var(sum(W5,2),1);

end
