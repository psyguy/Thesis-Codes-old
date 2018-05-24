clear global Network
global Network

disp('Preparing simulations')
tic 

% --- Time integration ---

Network.dt = 0.01; % Integration step (in ms)
simdur     = 2e3; % ms of simulation time
Network.nStep = round(simdur/Network.dt); % Number of time steps

Network.OutputVm       = 0; % 1=Output membrane potentials to Vm.txt, 0=no such output
Network.OutputSp       = 1; % 1=Output spike times to Spikes.txt, 0=no such output
Network.OutputSynState = 0; % 1=Output synapse state to SynState.txt, 0=no such output

fInit = fopen('Inits.txt','w');
fprintf(fInit,'%u\t%g\n%u\t%u\t%u',Network.nStep,Network.dt,Network.OutputVm,Network.OutputSp,Network.OutputSynState);
fclose(fInit);

% --- Populations ---
% generate Ncell cells, unconnected, that all have a constant synaptic
% input
Ncell = 400;
   Iinput = 0;
   gsyn = sort([linspace(-0.2,3,(Ncell-100)), linspace(-0.1,0.1,100)]);
AddPopulation(...
    'Name','Ex',...
    'nCell',Ncell,...
    'Type','HH',... 
    'Iinput',Iinput ,'Cm',10.0,'g_na',100.0,'g_naL',0.0175,'g_k',40.0,'g_kL',0.05,'g_clL',0.05,'phi',3.0,'E_k',Ek,'E_na',53,'E_cl',-82,...
    'gsyntest', gsyn, 'Esyntest', 50,...
    'V',-50,'n',0.07,'h',0.97,...
    'Position',repmat(linspace(0,10,Ncell)',1,3)...
);


% --- Output to files 
WriteCells;
WriteConnectivity;
toc

% --- run simulation in C++ environment
disp('Running simulations')
tic 
if(ismac)
    system('VerdandiLite_v1_0/VerdandiLite_OSX');
elseif(ispc)
    %system('VerdandiLite_v1_0\VerdandiLite_WIN32.exe');
    system('VerdandiLite_v1_0\VerdandiLite_CPP\BJZ.exe');
else
    error('No binaries available for VerdandiLite. Please compile C++ source files and add modify these Matlab lines.')
end
SkuldOut        % create output for analysis in Skuld (not used)
toc

Tana = [0,Inf]; % [ms] time span for analysis can be changed
Cellnums = (1:Ncell)-1;
method = 2; %use method 2, when input is constant, else use 1

Spikes = importdata('Spikes.txt');

% use only spikes in time span to analyze
ind = find(Spikes(:,1)>Tana(1) & Spikes(:,1)<Tana(2));
S = Spikes(ind,:);
S2 = sortrows(S,[2,1]);

switch method
    case 1      % average firing rate during entire simulation
        Nspike = hist(S(:,2),Cellnums);
        DT = (Tana(2)-Tana(1))/1e3;
        f = Nspike/DT;
    case 2      % spike rate at and of simulation ( 1/(last interspike interval) )
        f = zeros(size(Cellnums));
        Nspike = hist(S(:,2),Cellnums);
        for n = find(Nspike>2)
            lastspikepos = find(S2(:,2) == Cellnums(n),1,'last');            
            f(n) = 1e3/(S2(lastspikepos,1)-S2(lastspikepos-1,1)); % calculate difference in spike times between last and secondlast spike
        end
end
%% plot result
g = sort([linspace(-0.2,3,(Ncell*10-1000)), linspace(-0.1,0.1,1000)]);
F = interp1(gsyn,f,g,'pchip');
plot(gsyn,f,'r.'); hold on; xlabel('Synaptic conductane (mS/cm^2)'), ylabel('Spike rate')
plot(g,F); hold on; xlabel('Synaptic conductane (?/cm^2)'), ylabel('Spike rate')

%% save result
% 
