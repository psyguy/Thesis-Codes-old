clear all

% run detailed simulation using Brian, this script also contains the parameters of the simulation
% the simulation can take several minutes. To obtain the figure from the paper, set the simulated time to 2 seconds
system('brian_sim.py');

%% load simulation data from detailed model
load('gF_curve.mat','g','F')
load('Network_brian.mat')
Spikes(:,1) = Spikes(:,1)*1e3; % convert to ms
%% obtain and set parameters for simulation
tsim = t_brian;
Network.dt = diff(t_brian(1:2));

p.NoSTD = 0;    % 1 = set standard deviation of input to zero. To test effect of neglecting variability
p.F = [ones(1,100)*F(1),F,0];
p.g = [linspace(-.2,-0.0001,100),g,100];
p.std_cellparamEx = std(gsynex);
p.std_cellparamIn = std(gsynin);

% obtain input signal from network simulation
Fnoise = (hist(Spoisson,1:tsim(end)))/Ncellex*1e3;  % use ms wide bins
tnoise = (1:tsim(end))/1e3; %convert to seconds from ms
p.poisson = Fnoise;

p.g0 = g0; p.gamma = gamma; p.offsetg = offsetg;
p = parameters_meanfield(p,ConMat2,ConMat3,ConMat4,ConMat5);
p.T = [tsim(2), tsim(end)]/1e3; %convert to seconds from ms

%% run simulation
[sol,p] = meanfieldapprox(p);

%%
%% ------------------------------------------------
% code below is for plotting only
%%-------------------------------------------------

%% plot single cell firing rate curve
figure(320)
plot(g*(p.Eex-p.Vth),F)
xlabel('Input current (\muA/cm^2)')
ylabel('Firing rate (Hz)')

%% recalculate expressions (could be more efficient)
gsyntotEx = p.offsetg + sol.y(1,:) + sol.y(4,:)*p.Nsyn(2) + p.Cinh*sol.y(10,:)*p.Nsyn(4);
vargEx =  sol.y(3,:) + p.varN(2)*(sol.y(4,:)).^2 + p.Nsyna(2) * sol.y(6,:).^2 + p.varN(4)*(p.Cinh*sol.y(10,:)).^2 + p.Cinh^2*p.Nsyna(4) * sol.y(12,:).^2      ...
    + p.std_cellparamEx.^2;
if p.NoSTD; vargEx = 0*vargEx; end

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
gsyntotIn = p.offsetg + sol.y(7,:)*p.Nsyn(3) + p.Cinh*sol.y(13,:)*p.Nsyn(5);
vargIn =  p.varN(3)*(sol.y(7,:)).^2 + p.Nsyna(3) * sol.y(9,:).^2 + p.varN(5)*(p.Cinh*sol.y(13,:)).^2 + p.Cinh^2*p.Nsyna(5) * sol.y(15,:).^2      ...
    + p.std_cellparamIn.^2;   
    
if p.NoSTD; vargIn = 0*vargIn; end

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

%% compare synaptic conductances between NMM and network model, of all 5 synaptic connections
figure(321); subplot(2,2,1);
hold on
synInd = 1+[0, Ncellex,2*Ncellex, 2*Ncellex+Ncellin, 3*Ncellex+Ncellin, 3*Ncellex+2*Ncellin];
for count = 1
    plot(tsim/1e3,(p.Eex - p.Vth)*(mean(SynState(:,synInd(count):synInd(count+1)-1),2)),'cyan')
end
for count = 3
    plot(tsim/1e3,(p.Eex - p.Vth)*(mean(SynState(:,synInd(count):synInd(count+1)-1),2)),'green')
end
for count = 4
    plot(tsim/1e3,((p.Ein - p.Vth)*mean(SynState(:,synInd(count):synInd(count+1)-1),2)),'red')
end

%% std of synaptic conductances in network model, all 5 synaptic connections
subplot(2,2,2); hold on
count = 1;
plot(tsim/1e3,(p.Eex - p.Vth)*(std(SynState(:,synInd(count):synInd(count+1)-1),[],2)),'cyan')

for count = 3
    plot(tsim/1e3,(p.Eex - p.Vth)*(std(SynState(:,synInd(count):synInd(count+1)-1),[],2)),'green')
end
for count = 4
    plot(tsim/1e3,(abs((p.Ein - p.Vth))*std(SynState(:,synInd(count):synInd(count+1)-1),[],2)),'red')
end

%% firing rates of network model, both populations
subplot(2,2,3);hold on
Cellana= [0,Ncellex-1];
Ncellana = diff(Cellana)+1;
ind = Spikes(:,2)>=Cellana(1) & Spikes(:,2)<=Cellana(2);
ts = Spikes(ind,1)/1e3;
DTbin = 0.01;
[yhistEx,thistEx] = hist(ts,(DTbin/2:DTbin:tsim(end)/1e3));
f_networkEx = yhistEx/DTbin/Ncellana;
plot(thistEx,f_networkEx,'green.')

Cellana= [Ncellex,Ncellex+Ncellin-1];
Ncellana = diff(Cellana)+1;
ind = find(Spikes(:,2)>=Cellana(1) & Spikes(:,2)<=Cellana(2));
ts = Spikes(ind,1)/1e3;
DTbin = 0.01;
[yhistIn,thistIn] = hist(ts,(DTbin/2:DTbin:tsim(end)/1e3));
f_networkIn = yhistIn/DTbin/Ncellana;
plot(thistIn,f_networkIn,'red.')


%% plot NMM results

if p.NoSTD == false; plotcolor = 'blue'; else; plotcolor=[1 0 1]; end

%% synaptic conductances, of all 5 synaptic connections
subplot(2,2,1); hold on
plot(sol.x,(p.Eex - p.Vth)*sol.y(0*3+1,:),'--','Color',plotcolor)
%plot(sol.x,(p.Eex - p.Vth)*sol.y(1*3+1,:)*p.Nsyn(2),'--','Color',plotcolor)
plot(sol.x,(p.Eex - p.Vth)*sol.y(2*3+1,:)*p.Nsyn(3),'--','Color',plotcolor)
plot(sol.x,(p.Ein - p.Vth)*sol.y(3*3+1,:)*p.Nsyn(4),'--','Color',plotcolor)
%plot(sol.x,(p.Ein - p.Vth)*sol.y(4*3+1,:)*p.Nsyn(5),'--','Color',plotcolor)
legend('\eta->e','e->i','i->e')
xlabel('time (s)')
ylabel('I (\muA/cm^2)')
title('Mean synaptic current')
%% std of synaptic conductances, of all 5 synaptic connections
subplot(2,2,2); hold on
plot(sol.x,(p.Eex - p.Vth)*sqrt(sol.y(3,:)),'--','Color',plotcolor)
%plot(sol.x, (p.Eex - p.Vth)*sqrt(p.varN(2)*(sol.y(4,:)).^2 + p.Nsyna(2) * sol.y(6,:).^2),'--','Color',plotcolor )
plot(sol.x, (p.Eex - p.Vth)*sqrt(p.varN(3)*(sol.y(7,:)).^2 + p.Nsyna(3) * sol.y(9,:).^2),'--','Color',plotcolor )
plot(sol.x, abs((p.Ein - p.Vth))*sqrt(p.varN(4)*(sol.y(10,:)).^2 + p.Nsyna(4) * sol.y(12,:).^2),'--','Color',plotcolor )
%plot(sol.x, abs((p.Ein - p.Vth))*sqrt(p.varN(5)*(sol.y(13,:)).^2 + p.Nsyna(5) * sol.y(15,:).^2),'--','Color',plotcolor )
legend('\eta->e','e->i','i->e','NMM')
xlabel('time (s)')
ylabel('\sigma_I (\muA/cm^2)')
title('Std of synaptic currents')

%% firing rates, of both populations
subplot(2,2,3);hold on
plot(sol.x,fEx,'--','Color',plotcolor);
plot(sol.x,fIn,'--','Color',plotcolor);
xlabel('time (s)')
ylabel('f (Hz)')
title('Mean firing rates')
legend('f_e','f_i','NMM')

%% std of firing rates, of both populations
subplot(2,2,4); hold on

Cellana= [0,Ncellex-1];
F_inst_sim_ex = zeros(length(sol.x),Ncellana);
for c1 = Cellana(1):Cellana(2);
    startpos = 1;
    xposold = 1;
    for c2 = find(Spikes(:,2)==c1)';
        
        [~,xposnew] = min(abs(sol.x*1e3-Spikes(c2,1)));
                
        F_inst_sim_ex(xposold:xposnew,c1-Cellana(1)+1) = 1e3./(Spikes(c2,1)-Spikes(startpos,1));
        startpos = c2+1;
        xposold = xposnew;
    end
end
plot(sol.x, std(F_inst_sim_ex,0,2),'g'); 

Cellana= [Ncellex,Ncellex+Ncellin-1];
F_inst_sim_in = zeros(length(sol.x),Ncellana);
for c1 = Cellana(1):Cellana(2);
    startpos = 1;
    xposold = 1;
    for c2 = find(Spikes(:,2)==c1)';
        [~,xposnew] = min(abs(sol.x*1e3-Spikes(c2,1)));
        F_inst_sim_in(xposold:xposnew,c1-Cellana(1)+1) = 1e3./(Spikes(c2,1)-Spikes(startpos,1));
        startpos = c2+1;
        xposold = xposnew;
    end
end
plot(sol.x, std(F_inst_sim_in,0,2),'r'); 

plot(sol.x,sqrt(varfEx),'--','Color',plotcolor);
plot(sol.x,sqrt(varfIn),'--','Color',plotcolor);
xlabel('time (s)')
ylabel('\sigma_f (Hz)')
title('Std of firing rates')
legend('f_e','f_i','NMM')
