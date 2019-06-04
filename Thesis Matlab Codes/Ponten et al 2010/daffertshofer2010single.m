function [E, t] = daffertshofer2010single(varargin)
% p_ext, n, t0, E_k, E_prev
% Make sure if you only want to connect them together you need to
% explicitly assign n=1 and E_k in the function arguments

%% defining parameters
A_e = 1.6; % for exitatory 
a_e = 55;
b_e = 605;
A_i = 32; % for inhibitory
a_i = 27.5;
b_i = 55;

G1 = 32;
G2 = 3;

% pMu = 550;
% pSig = 0.1; % is it correct?
% p_ext = normrnd(pMu,pSig);

%% Init

defaults = {0,1,1,rand(1),rand(1)};
defaults(1:4) = varargin;

p_ext = cell2mat(defaults(1));
n = cell2mat(defaults(2));
t0 = cell2mat(defaults(3));
E_k = cell2mat(defaults(4));
E_prev = cell2mat(defaults(5));



% we first calculate the inhibitory population since it makes I_j for us
% hence we do not need to pass it to through the argument everytime we call
% the function

% inhibitory interneurons, the bottom half of the model
    V_i_prev = G1*E_prev*hfunc(t0,A_e,a_e,b_e);
    I_prev = sfunc(V_i_prev);
    
for (t=t0:(t0+n))
    
    E = E_prev;
    I = I_prev;
    % exitatory (pyramidal) neurons
    V_en = (E_k+p_ext)*hfunc(t,A_e,a_e,b_e) + I*G2*hfunc(t,A_i,a_i,b_i);
    E_n = sfunc(V_en);
    
    % inhibitory interneurons, the bottom half of the model
    V_in = G1*E*hfunc(t,A_e,a_e,b_e);
    I_jn = sfunc(V_in);

    E(t+1) = E_n;
    
    E_prev = E_n
    I_prev = I_jn

end

end