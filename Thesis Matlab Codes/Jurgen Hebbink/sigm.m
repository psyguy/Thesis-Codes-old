%Computes Sigmoidal firing function in the Jansen-Rit blocks of the
%model. Parameter values are hard coded.
function outp=sigm(v)
    outp=5./(1+exp(0.56*(6-v)));
end

