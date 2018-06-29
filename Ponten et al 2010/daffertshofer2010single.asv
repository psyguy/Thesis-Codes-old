function [E_j,I_j] = daffertshofer2010single(t,E_k,p_ext,E_j_prev)
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

% we first calculate the inhibitory population since it makes I_j for us
% hence we do not need to pass it to through the argument everytime we call
% the function
E_j = E_j_prev;
% inhibitory interneurons
V_i = G1*E_j*hfunc(t,A_e,a_e,b_e);
I_j = sfunc(V_i);

% exitatory (pyramidal) neurons
V_e = (E_k+p_ext)*hfunc(t,A_e,a_e,b_e) + I_j*G2*hfunc(t,A_i,a_i,b_i);
E_j = sfunc(V_e);

end

%%
function h = hfunc(t,A,a,b)

if (t<0)
    h = 0;
else
    h = A*(exp(-a*t)-exp(-b*t));
end

end

%%
function sv = sfunc(v)

q = 0.34;
v0 = 7;
g = 25;

if (v<v0)
    sv = g*exp(q*(v-v0));
else
    sv = 2*g - g*exp(-q*(v-v0));
end

end


