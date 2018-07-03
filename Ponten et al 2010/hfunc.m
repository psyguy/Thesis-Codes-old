function h = hfunc(t,A,a,b)

if (t<0)
    h = 0;
else
    h = A*(exp(-a*t)-exp(-b*t));
end

end
