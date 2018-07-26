%-------------------------------------------------------------------------%
%-                           Spase Petkoski                              -%
%- Theoretical Neuroscience Group -- Institute of Systemes Neuroscience  -%
%-                        2018, AMU Marseille                            -%
%-------------------------------------------------------------------------%

%------- connectome  with N x tau_m values in the phases
%-------
% Script for calculating phase dynamics for one subject and PLV related
% metrics and statistics.
% Each node has a same natural frequency (fm) and coupling strenghts and 
% link lenghts are from the connectome. Time delays are obtained assuming
% constant propagation velocity (v_cnctm). 

clc
clear all
close all
noiseseed=1;
rng(noiseseed)
global mu fm K_cnctm h dt N Nt sc Nds t D K K0 v_cnctm tfin l_cnctm T prntfgs
hpar=0.05; % parameter for the time step; guarantees that h*max(f(x)) < hpar
%----control parameters-----
incoh=2; % 0 for all equal; 1 for complete incoherence; (0,1) partial incoherence; (1,2] same as (0,1], but the hemispheres are in antiphase
D=2; % noise intensity
K0=0.05; % fx0l1w2f40t120K3D20v6sc1h005; fx0l1w2f20t120K05D1v5sc5h005.png
scalweigh=0; % 1 for log, 0 for no scaling; 2 for sqrt, 4 for ^1/4
tfin0=120; % total time of the simulation (in seconds)
fm=8; % natural frequency of the oscillators in Hz
v_cnctm=5000; % mm/s
prntfgs=1; % to print figures or no
%%
Nfmax=2; % minimum 2; determines the step dt; dt=1/fmax=1/(Nfmax*fnyquist)
WLp=10; % number of slowest cycles in the windows for plv
overlap=0.75; %0.5; % overlapping of windows
cnt0=20; % counter at every cnt seconds
%
%-------- connectome parameters ---------------
load('wl68N100.mat') % the file contains connectomes of 100 subjects from HCP
subj=50; % which subject
K_cnctm_raw=squeeze(wght(:,:,subj)); % strengths of the links 
l_cnctm=squeeze(lnght(:,:,subj)); % lengths of the links
N=size(K_cnctm_raw,1);

%% sorting links/nodes by strength
% sort the inter and itrahemisphere links by strength
[Bind,Kind]=sort(sum(K_cnctm_raw));
Kindl=Kind(find(Kind<35));
Kindr=Kind(find(Kind>34));
Ikin=[Kindl, Kindr]; % nodes sorted by instrength
%%
tfin=tfin0;
T=round(min(max(round(tfin/4), 20),50)); % initial transition time to be discarded.
mu=fm*2*pi;  % natural frequency of the oscillators in rad/s
%
K_cnctm_raw(logical(eye(size(K_cnctm_raw))))=0;
disp(['f= ', num2str(fm) 'Hz', ', K0= ', num2str(K0), ', v= ', num2str(v_cnctm), ', D= ', num2str(D), ', K0= ', num2str(K0), ', hpar= ', num2str(hpar)])
%
sc=length(K_cnctm_raw(K_cnctm_raw>0))/N; % scaling parameter for the sum of the couplings; average connectivity
Km1=median(K_cnctm_raw(:));
K=1/Km1/sc; % scaling taken into account in the coupling strenght; for all-to-all this would correspond to K=K/N
K_cnctm=K0*K*K_cnctm_raw;
%--- setting the time step ------
h=hpar/(max([max(K_cnctm(:))/10, mu*5*hpar, D/100, 1])); % at least 5 points per cycle. , max(sum(K_cnctm))
dt=floor(pi/Nfmax/mu/h)*h ; % fNyqust=2*fmax=Nfmax*mu/pi; dt=pi/Nfmax/mu~3/Nfmax/mu; for N=3.14 times larger frequency than f, 1 is obtained
if dt<h, dt=h; end;
disp(['h=', num2str(h), '; dt=', num2str(dt), '; Kmax=', num2str(max(K_cnctm(:)), '%2.2f'), ...
    '; Kmean=', num2str(mean(K_cnctm(find(triu(ones(N),1)))), '%3.3f'), '; v_ef=', num2str(round(v_cnctm))])
t=0:dt:tfin+2*dt; % so that it is ensured that tmax>=tfin; otherwise, if tfin+dt, with floor the last rounding (downsampling) might finish before tfin
Nt=length(t);
Nh=floor((tfin+2*dt)/h); %
Nds=dt/h; % downsampling points, % there is a problem if h is not a multiple of dt
if fix(Nds)~=Nds && abs(round(Nds)-Nds>1e-12),
    error('h should be a multiple of dt.');
else
    Nds=round(Nds);
end;
%---
l_cnctm(K_cnctm_raw==0)=0;
tau_cnctm=round(l_cnctm/v_cnctm/h); % matrix with time delays in units of time-steps of the integration
tau_m=max(max(tau_cnctm(:)), Nds)+1; % maximum time-delay used for setting the length of the buffer
tau1 = round(mean(wtl([tau_cnctm(1:N/2,1:N/2), tau_cnctm(1+N/2:N,1+N/2:N)], [K_cnctm_raw(1:N/2,1:N/2), K_cnctm_raw(1+N/2:N,1+N/2:N)],0))); % mean intrahemispheric delay
tau2 = round(mean(wtl(tau_cnctm(1:N/2,1+N/2:N), K_cnctm_raw(1:N/2,1+N/2:N),0))); % mean interhemispheric delay
%
clear K_cnctm_raw;
disp(['tau_1= ', num2str(1000*tau1*h),'ms; tau_2= ', num2str(1000*tau2*h),'ms; T= ', num2str(1000/fm), 'ms'])
%% allocation and obtaining indices of delayed phases
[theta, thds, i]=alloc_cnctmPlosCB(incoh, tau_m);
[idx, ~, ~, ~, ~, idx1]=ind_cnctmPlosCB(theta, tau_cnctm);
%%
varstrKM.idx=idx;
varstrKM.idx1=idx1;
varstrKM.mu=mu;
varstrKM.K_cnctm=K_cnctm;
varstrKM.h=h;
varstrKM.Nds=Nds;
varstrKM.t=t;
varstrKM.cnt=round(cnt0/dt);
varstrKM.N=N;
varstrKM.N2=N/2;
%----
cnt=round(cnt0/dt);
eta=sqrt(2*h*D)*randn(N, Nh-i); % i=tau_m
while i<Nh
    %%%% Heun integration scheme %%%%
    chck=ceil((i+1)/Nds);   % counter for the downsampled phases until chck=(i+1)/Nds, all the thds and thds2 are returned NaN
    [theta,  thds(:,chck+1),  i, ~] = KMcnctmHfPlosCB(theta,   eta(:,i+1-tau_m), i,   chck, varstrKM);
end
% still to be able to print in case of breaking the simulation
[~,I]=find(isnan(thds), 1, 'first');
if isempty(I),
    I=length(thds(1,:))+1;
else
    if I<length(t)-1
        warning('NaN detected|')
        disp(['I=', int2str(I), ',  t=', num2str(t(I))])
    end
end;
t(I:end)=[];
Nt=length(t);
tfin=floor(t(I-1));
thds(:,I:end,:)=[];
%% -------------
% printing
fname=['PlosCBsubj' int2str(subj) 'f' int2str(fm)  't'  int2str(tfin)...
    'K' num2str(K0) 'D' num2str(D) 'v' num2str(v_cnctm/1000) 'h' num2str(hpar)];
fname=regexprep(fname,{'\.','txt'},{'','.txt'});
%[zs, zs1, zs2]=printFTevol(thds, fname, 0, 0.05); % print time-series of order parameters
%
zs = mean(exp(1j*thds),1); %  complex order parameters (whole brain)
zs1 = mean(exp(1j*thds(1:N/2,:))); %  complex order parameters (left hemisphere)
zs2 = mean(exp(1j*thds(N/2+1:N,:))); %  complex order parameters (right hemisphere)
Omega = diff(unwrap(angle(zs)))/dt/2/pi;  % frequency of synchronization for each of the hemispheres (time series)
Omega_mean=mean(Omega(round(Nt/2):end)); % mean frequency of synchronization for each of the hemispheres (makes sense for stationary solutions)
%% variables to be plot for the entrained oscillators
ki=sum(K_cnctm); % in-strength of the nodes
taul=[tau1, tau2]*h; % used for calculating \Omega \Delta \tau and \tilde{\Omega} \tau
% printing
relativeph=1; % if 1 then the phases in the right hemisphere are plot relative to those in the left
printsep=2; % print figures separately
%printKi_ph(thds, zs1(:), zs2(:), mu, D, ki, taul, dt, t, K0, fname, relativeph, printsep)
%
%% PLV, phase lags
% -- calculate PLV and statistics
pi2=2*pi;
WL=WLp./min([mean(diff(unwrap(angle(zs1(:,1))))/dt/2/pi), mean(diff(unwrap(angle(zs2(:,1))))/dt/2/pi),fm]);
%
[R,xR,surr,surr2,thsync2,thsyncplv2,thsync,thsyncplv,Dthmean,Dthvar,Dplvmean,...
    Dplvvar,Dthmean2,Dthvar2,Dplvmean2,Dplvvar2,dth, PLVs, phiplv]= ...
    plvstatPlosCB(thds, WL, overlap, []);
Rm=squeeze(mean(R,3));

% %%
% % print statisitcs of PLV metrics
% fname=[fname 'O' int2str(round(Omega_mean)) 'w' int2str(WLp) 'k' int2str(overlap)];
% %
% printDph(fname, surr, surr2, Dthmean,Dthvar,Dplvmean,Dplvvar,Dthmean2,Dthvar2,Dplvmean2,Dplvvar2, ...
%     [],[],[],[],[],[],[],[], NaN, R, [], PLVs, []);
% 
% %% -- print pairwise time-series --
% xlim1=round(tfin/2) - 5; %0.05*round(tfin);
% xlim2=round(tfin/2) + 5; %0.05*round(tfin);
% %% pairs by the strength of the connection
% m=5; % nodes for which to be printed
% n=10; % nodes for which to be printed
% printDph_ch(fname, m, n, tau1, tau2, 1, squeeze(dth(m,n,:)), surr, surr2(m,n),  PLVs(m,n), squeeze(R(m,n,:)),...
%     squeeze(phiplv(m,n,:)), xR, squeeze(thsync(m,n,:)), squeeze(thsync2(m,n,:)), squeeze(thsyncplv(m,n,:)),...
%     squeeze(thsyncplv2(m,n,:)), 0, [int2str(m) '/' int2str(m)], xlim1, xlim2);
% %----

