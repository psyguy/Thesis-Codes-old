%Computes Sigmoidal firing function in the Jansen-Rit blocks of the
%model. Parameter values are hard coded.

% Edited by Manuel 20180629@1314
function outp=sigm(v,pars)
    v0 = pars(1);
    vmax = pars(2);
    r = pars(3);
    outp=vmax./(1+exp(r*(v0-v)));
end

