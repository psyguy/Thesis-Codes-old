function [E_j] = daffertshofer2010single(varargin)
% p_ext, n, E_k, E_j_prev
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

defaults = {0,1,rand(1),rand(1)};
defaults(1:nargin) = varargin;

p_ext = cell2mat(defaults(1));
n = cell2mat(defaults(2));
E_k = cell2mat(defaults(3));
E_j_prev = cell2mat(defaults(4));



% we first calculate the inhibitory population since it makes I_j for us
% hence we do not need to pass it to through the argument everytime we call
% the function

for (t=1:n)
    
    E_j = E_j_prev;
    % inhibitory interneurons, the bottom half of the model
    V_in = G1*E_j*hfunc(t,A_e,a_e,b_e);
    I_jn = sfunc(V_in);

    % exitatory (pyramidal) neurons
    V_en = (E_k+p_ext)*hfunc(t,A_e,a_e,b_e) + I_jn*G2*hfunc(t,A_i,a_i,b_i);
    E_jn = sfunc(V_en);
    
    E_j(t) = E_jn;
    E_j_prev = E_jn;

end

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


