function [SolWC,p] = meanfieldapprox(p)
%warning('no std taken into account');
if p.solvewithode23s
    p.reltoll = 1e-3;
    p.maxstep = 1;
    tspan = p.T; options = odeset('RelTol',p.reltoll,'MaxStep',p.maxstep);%,'OutputFcn',@odeprog,'Events',@odeabort );
    funhan = (@(t,h)WC(t,h,p));
    SolWC = ode23s(funhan,tspan,[p.h0],options); %  Solve the differential equations
else
    % t in seconds
    t = p.T(1):p.Dt:(p.T(2)-p.Dt);
    Nt = length(t); %ceil((p.T(2)-p.T(1))/p.Dt);
    h = zeros(length(p.h0), Nt);
    h(:,1) = p.h0;
    
    for nt = 1:Nt-1
        dh = WC(t(nt),h(:,nt),p);
        h(:,nt+1) = h(:,nt) + dh*p.Dt;% + sqrtdt*noise;
        if ~mod(nt,round(Nt/30))
            fprintf('.')
        end
    end
    
    SolWC.y = h;
    SolWC.x = t;
end
end

% Derivatives
function dh = WC(t,h,p)  %Wilson-Cowan type model
dh = zeros(15,1);
e1 = exp(1);

%npoisson = ceil(t*1e3);
lambda_ext = p.lambda; %can be rewritten to put in a time dependent signal

% Calculate variances in total synaptic input and frequency output
% excitatory population
gsyntotEx = h(1)           +  h(4)*p.Nsyn(2)                                    +  p.Cinh * h(10) *p.Nsyn(4);      % synapses 1,2 and 4
vargEx    = h(3)           + p.varN(2)*(h(4)).^2 + p.Nsyna(2) * h(6).^2  + p.varN(4)*(p.Cinh*h(10)).^2 + p.Cinh^2*p.Nsyna(4) * h(12).^2      ...
    + p.std_cellparamEx.^2;  % variance due to poisson input, variance in population input,  and different number of synapses respectively

mu = gsyntotEx;
width = vargEx;
if width < 1e-10;   % prevent devision by zero errors.
    f(1) = interp1(p.g,p.F,gsyntotEx,'linear','extrap');
    varf(1) = 0;
else
    % f(mu,width) and varf(mu,width) could maybe be tabulated to improve efficiency
    gauss = exp(-(p.g-mu).^2/(2*width))/sqrt(2*pi*width);
    f(1) = trapz(p.g,gauss.*p.F);
    varf(1) = trapz(p.g,gauss.*(p.F-f(1)).^2);
end
% inhibitatory population
gsyntotIn  = h(7)*p.Nsyn(3)                             + p.Cinh*h(13) *p.Nsyn(5); % synapses 3 and 5
varg_In    = p.varN(3)*(h(7)).^2 + p.Nsyna(3) * h(9).^2 + p.varN(5)*(p.Cinh*h(13)).^2 + p.Cinh^2*p.Nsyna(5) * h(15).^2  ...
    + p.std_cellparamEx.^2;        % variance due to variance in population input, noisy input and different number of synapses respectively

mu = gsyntotIn;
width = varg_In;   % prevent devision by zero errors.
if width < 1e-10;
    f(2) = interp1(p.g,p.F,gsyntotEx,'linear','extrap');
    varf(2) = 0;
else
    gauss = exp(-(p.g-mu).^2/(2*width))/sqrt(2*pi*width);
    f(2) = trapz(p.g,gauss.*p.F);
    varf(2) = trapz(p.g,gauss.*(p.F-f(2)).^2);
end

% Synaptic dynamics
%h(3n+1) = g
%h(3n+2) = dg/dt
%h(3n+3) = sigma_g

% poisson firing synapse (extracortical input -> E cells)
count = 1;
n = (count-1)*3+1;
dh(n) = h(n+1);              
dh(n+1) = -2*p.gamma(count).*h(n+1)  - p.gamma(count).^2.*h(n) + e1*p.g0(count).*p.gamma(count).*lambda_ext;
dh(n+2) = -2*h(n+2)/p.tau_approx(1) + p.g0_approx(1)^2*lambda_ext; % lambda convoluted with H^2

% regular firing synapses
origin = [NaN, 1, 1, 2, 2]; % postsynaptic population, 1 = e, 2 = i
for count = 2:5
    n = (count-1)*3+1;
    dh(n) = h(n+1);         
    dh(n+1) = -2*p.gamma(count).*h(n+1)  - p.gamma(count).^2.*h(n) + e1*p.g0(count).*p.gamma(count).*f(origin(count));
    dh(n+2) = -h(n+2)/(p.tau_approx(count)) + p.g0_approx(count)*sqrt(varf(origin(count))); %sigma_f convoluted with H
end

end