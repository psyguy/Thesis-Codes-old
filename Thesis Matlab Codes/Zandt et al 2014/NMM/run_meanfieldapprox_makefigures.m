clear all

%% load F(g)-curve
% the F(g)-curve is used, rather than the F(I)-curve.
% This is slightly mor accurate.
load('gF_curve.mat','g','F')
p.F = [ones(1,100)*F(1),F,0];
p.g = [linspace(-.2,-0.0001,100),g,100];

%% load parameters
p = parameters_meanfield(p);

% simulation parameters are stored in structure p:
% F(g), g0(1:5), gamma(1:5), Cinh,
% std_cellparamEx, std_cellparamIn
% Nsyn(2:5), Nsyna(2:5), k(2:5)
% tau_approx(2:5), g0_approx(2:5)                % time constants and amplitudes of exponential approximations

%% Run actual simulation
[sol,p] = meanfieldapprox(p);


%% ------------------------------------------------
% code below is for plotting only
%%-------------------------------------------------

%% recalculate expressions (could be more efficient)

gsyntotEx = sol.y(1,:) + sol.y(4,:)*p.Nsyn(2) + p.Cinh*sol.y(10,:)*p.Nsyn(4);
vargEx =  sol.y(3,:) + p.varN(2)*(sol.y(4,:)).^2 + p.Nsyna(2) * sol.y(6,:).^2  + p.varN(4)*(p.Cinh*sol.y(10,:)).^2 + p.Cinh^2*p.Nsyna(4) * sol.y(12,:).^2      ...
    + p.std_cellparamEx.^2;   
    
m = gsyntotEx;
fEx = zeros(1,length(m));
varfEx = zeros(size(fEx));
for count = 1:length(m)
if vargEx(count) < 1e-10;   % prevent devision by zero errors.
    fEx(count) = interp1(p.g,p.F,gsyntotEx(count),'linear','extrap');
    varfEx(count) = 0;
else
    gauss = exp(-(p.g-m(count)).^2/(2*vargEx(count)))/sqrt(2*pi*vargEx(count));
    fEx(count) = trapz(p.g,gauss.*p.F);
    varfEx(count) = trapz(p.g,gauss.*(p.F-fEx(count)).^2);
end
end
gsyntotIn = sol.y(7,:)*p.Nsyn(3) + p.Cinh*sol.y(13,:)*p.Nsyn(5);
vargIn =  p.varN(3)*(sol.y(7,:)).^2 + p.Nsyna(3) * sol.y(9,:).^2  + p.varN(5)*(p.Cinh*sol.y(13,:)).^2 + p.Cinh^2*p.Nsyna(5) * sol.y(15,:).^2      ...
    + p.std_cellparamIn.^2;   
    
m = gsyntotIn;
fIn = zeros(1,length(m));
varfIn = zeros(size(fIn));
for count = 1:length(m)
if vargIn(count) < 1e-10;
    fIn(count) = interp1(p.g,p.F,gsyntotIn(count),'linear','extrap');
    varfIn(count) = 0;
else
    gauss = exp(-(p.g-m(count)).^2/(2*vargIn(count)))/sqrt(2*pi*vargIn(count));
    fIn(count) = trapz(p.g,gauss.*p.F);
    varfIn(count) = trapz(p.g,gauss.*(p.F-fIn(count)).^2);
end
end

%% plot single cell firing rate curve
figure(322)
plot(g*(p.Eex-p.Vth),F)
xlabel('Input current (\muA/cm^2)')
ylabel('Firing rate (Hz)')


%% synaptic conductances, of all 5 synaptic connections
figure(323);
subplot(2,2,1); hold on
plot(sol.x,(p.Eex-p.Vth)*sol.y(0*3+1,:),'blue')
plot(sol.x,(p.Eex-p.Vth)*sol.y(1*3+1,:)*p.Nsyn(2),'green')
plot(sol.x,(p.Eex-p.Vth)*sol.y(2*3+1,:)*p.Nsyn(3),'green--')
plot(sol.x,(p.Ein-p.Vth)*sol.y(3*3+1,:)*p.Nsyn(4),'red--')
plot(sol.x,(p.Ein-p.Vth)*sol.y(4*3+1,:)*p.Nsyn(5),'red')
xlabel('time (s)')
ylabel('I (\muA/cm^2)')
title('Mean synaptic current')
legend('external \eta->e','e->e','e->i','i->e','i->i')
xlabel('time (s)')
%% std of synaptic conductances, of all 5 synaptic connections
subplot(2,2,2); hold on
plot(sol.x, (p.Eex-p.Vth)*sqrt(sol.y(3,:)),'blue')
plot(sol.x, (p.Eex-p.Vth)*sqrt(p.varN(2)*(sol.y(4,:)).^2 + p.Nsyna(2) * sol.y(6,:).^2),'green' )
plot(sol.x, (p.Eex-p.Vth)*sqrt(p.varN(3)*(sol.y(7,:)).^2 + p.Nsyna(3) * sol.y(9,:).^2),'green--' )
plot(sol.x, abs(p.Ein-p.Vth)*sqrt(p.varN(4)*(sol.y(10,:)).^2 + p.Nsyna(4) * sol.y(12,:).^2),'red--' )
plot(sol.x, abs(p.Ein-p.Vth)*sqrt(p.varN(5)*(sol.y(13,:)).^2 + p.Nsyna(5) * sol.y(15,:).^2),'red' )
legend('external \eta->e','e->e','e->i','i->e','i->i')
xlabel('time (s)')
ylabel('\sigma_I (\muA/cm^2)')
title('Std of synaptic currents')

%% mean population firing rates
subplot(2,2,3);hold on
plot(sol.x,fEx,'green');
plot(sol.x,fIn,'r');
title('Mean firing rates')
xlabel('time (s)')
ylabel('firing rate (Hz)')
legend('f_e','f_i')
%% std of firing rates in population
subplot(2,2,4);hold on
plot(sol.x,sqrt(varfEx),'green');
plot(sol.x,sqrt(varfIn),'r');
xlabel('time (s)')
ylabel('std of firing rate (Hz)')
title('std of firing rates in population')
legend('\sigma_{f,e}','\sigma_{f,i}')