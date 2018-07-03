function sv = sfunc(v)

q = 0.34e3;
v0 = 7e-3;
g = 25;

if (v<v0)
    sv = g*exp(q*(v-v0));
else
    sv = 2*g - g*exp(-q*(v-v0));
end

end