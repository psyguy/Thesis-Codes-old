options = optimoptions('fsolve','Display','none');%,'PlotFcn',@optimplotfirstorderopt);

fun = @grifau2006simple;
y0 = [0,0,0,0];
y = fsolve(fun,y0,options)


%% Trying manually
p = 120;
yinit = [0,0,0,0];
Y(1) = yinit;
for (i=1:100)
    Y(i+1,:) = grifau2006simple(p, Y(i,:));
%     Yprev = Y;
end
