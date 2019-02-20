%-------------------------------------------------------------------------%
%-                           Spase Petkoski                              -%
%- Theoretical Neuroscience Group -- Institute of Systemes Neuroscience  -%
%-                        2018, AMU Marseille                            -%
%-------------------------------------------------------------------------%

%------- network with N nodes and buffer of length tau_m for the phases
% script for calculating phase dynamics and statistics for N oscillators
% connected all-to-all (although zero weights are also allowed)
%
% freqnoise: 1|2
%                - 1 distributed natural frequencies (Lorentzian).  
%                   mu and gamma are the mean and the width
%                - 2 noise (white Gaussian)
%                   D is the standard deviation of the noise
% Klog: 0|1
%                - 0 homogeneous coupling weights K0
%                - 1 log distributed weights with mean K0 and variance 2*K0
%                (can be modified)
% NA: 0|1
%                - 0 fixed parameters (natural frequencies and couplings) 
%                - 1 non-autonomous (NA) natural frequencies, only for N=2
%                and the frequencies are equal (doesn't make sense to have
%                them distributed when N=2) with or without noise. Coupling
%                strengths are also modulated in this case. The natural 
%                frequencies are modulated with strength etaO and frequency 
%                OmegaO, while for the couplings they are Omega and eta.
%                effective frequencies Oeff and couplings Keff appear 
% Arnd: 1|2 
%                - 1 spatialy ordered time delays 
%                - 2 random (spatialy homogeneous) time-delays                
%
% omega         natural frequencies at each node

clc
clear all
close all
global mu fm K_cnctm h dt N Nt omega  Nds t D K K0 tfin l_cnctm par T 
rng default; 
hpar=0.02; % parameter for the time step; guarantees that h*max(f(x)) < hpar
incoh=0; % initial coherence of the phases, 0 for incoherent (equaly distirbuted) 
Nfmax=2; % minimum 2; determines the step dt; dt=1/fmax=1/(Nfmax*fnyquist)
%----control parameters-----
D=5; % noise intensity; it corresponds to gamma for deterministic case
K0=30; % coupling strength
Klog=0; % log-normal K or no
tfin0=200; % total time of the simulation (in seconds)
fm=4; % mean of the natural frequencies in Hz
taul=[0.03, 0.65]; % delays in seconds
N=200; % number of oscillators
NA=1; % 0 for no non-autonomicity and 1 for non-autonomous omega and/or K
if N==2, taul =taul(1); else NA=0; end
domega=0.03; % relative difference between the means of each cluster, used only for N=2
%%
%
freqnoise=1; % 1 for noise, 0 for frequency heterogeneity
Arnd=2; % 1 for ordered time delays, 2 for random
if freqnoise, Klog=1; end % coupling are fixed in all cases of distributed frequencies 
%
cnt0=100; % counter at every cnt seconds
%-------- connectome parameters ---------------
K_cnctm=K0*ones(N);
if N>2
    Isupdiag    = find(triu(ones(N),+1)); % Indexes of all the values above the diagonal.
    if Klog
        mK = K0; % mean
        vK = 2*K0; % variance
        muK = log((mK^2)/sqrt(vK+mK^2));
        sigmaK = sqrt(log(vK/(mK^2)+1));
        [M,V]= lognstat(muK,sigmaK);
        
        Kl = lognrnd(muK,sigmaK,N/2,1);
        Kr=Kl(randperm(N/2));
        Kl=Kl(randperm(N/2));
        Ki=[Kl;Kr];
        for i=1:N
            rnd=rand(N,1); % uniformly distributed
            K_cnctm(i,randperm(N))= Ki(i) + 0*Ki(i)*2*(rnd - mean(rnd)); % adding rand part doesn;t seem to change anything in Ki so better not to; however makes the distribution of K_cnctm wider
        end
        K_cnctm=(K_cnctm+K_cnctm')/2; 
    end
end
K_cnctm=K_cnctm/N;
K_cnctm(logical(eye(size(K_cnctm))))=0;
%
if ~freqnoise
    gamma=D;
    D=0;
    mu=2*pi*fm;
    if N>2
        omega1(1:N/2,1) = mu + tan(linspace(-pi/2+0.000001, pi/2-0.000001, N/2)').*gamma;
        omega2(1:N/2,1) = mu + tan(linspace(-pi/2+0.000001, pi/2-0.000001, N/2)').*gamma;
    else
        omega1 = mu; omega2 = mu; % might be overwritten below
    end
    omega=([omega1; omega2]);
    rp=randperm(N);
end
% create ordered and random 
if N>2
    % ordered delays
    C11 = tril(ones(N/2),-1)+triu(ones(N/2),1);
    C44 = tril(ones(N/2),-1)+triu(ones(N/2),1);
    C22 = ones(N/2,N/2);
    C33 = C22';
    %-----
    C1=taul(1)*C11;
    C2=taul(2)*C22;
    C3=taul(2)*C33;
    C4=taul(1)*C44;
    %----
    tau5ord= [[C1, C2]; [C3, C4]]; % realistic, ordered as tau1 tau2 tau2 tau1
    % random
    tau5 = triu(ones(N),1);
    p1=0.5;
    Nnd=N*(N-1)/2; % number of links in the uper triangle
    Nnd1=round(p1*Nnd); % number of links with tau1 in the upper triangle
    tau_nd(1:Nnd1)=taul(1);
    tau_nd(Nnd1+1:Nnd)=taul(2);
    tau_nd(randperm(Nnd))=tau_nd;
    tau5(tau5==1)=tau_nd;
    tau5=tau5+tau5';
    if Arnd==1
        l_cnctm=tau5ord;
        tit='A';
    elseif Arnd==2
        l_cnctm=tau5;
        tit='rnd';
    end
else
    l_cnctm=[0,taul;taul,0];
    tit='2 osc';
end
%
disp(['mod ', tit, ', D= ', num2str(D), ', K0= ', num2str(K0), ', hpar= ', num2str(hpar)])
if N>2
    ki=sum(K_cnctm,2); 
end
%
tfin=tfin0;
T=round(min(max(round(tfin/4), 20),50)); % initial transition time
%
mu=fm*2*pi;
disp([ 'f= ', num2str(fm) 'Hz', ', K0= ', num2str(K0)])
%%
%--- setting the time step ------
h=hpar/(max([max(K_cnctm(:)), mu*5*hpar, D, 1])); % at least 5 points per cycle. , max(sum(K_cnctm))
dt=floor(pi/Nfmax/mu/h)*h; % fNyqust=2*fmax=Nfmax*mu/pi; dt=pi/Nfmax/mu~3/Nfmax/mu; for N=3.14 times larger frequency than f, 1 is obtained
if dt<h, dt=h; end;
disp(['h=', num2str(h), '; dt=', num2str(dt), '; Kmax=', num2str(max(K_cnctm(:)), '%2.2f'), ...
    '; Kmean=', num2str(mean(K_cnctm(find(triu(ones(N),1)))), '%3.3f')])
t=0:dt:tfin+2*dt; % so that it is ensured that tmax>=tfin; otherwise, if tfin+dt, with floor the last rounding (downsampling) might finish before tfin
Nt=length(t);
Nh=floor((tfin+2*dt)/h); %
Nds=dt/h; % downsampling points, % there is a problem if h is not a multiple of dt
if fix(Nds)~=Nds && abs(round(Nds)-Nds>1e-12),
    error('h should be a multiple of dt.');
else
    Nds=round(Nds);
end
l_cnctm(K_cnctm==0)=0;
tau_cnctm=round(l_cnctm/h);
tau_m=max(max(tau_cnctm(:)), Nds)+1;
clear K_cnctm_raw_round;
%% allocation and obtaining indices of delayed phases
[theta, thds, i]=alloc_cnctmPlosCB(incoh, tau_m); % also returns i and states when simulations starts
[idx, ~, ~, ~, ~, idx1]=ind_cnctmPlosCB(theta, tau_cnctm);
%%
varstrKM.K_cnctm=K_cnctm;
if N>2
    if freqnoise
        varstrKM.omega=mu;
    else
        varstrKM.omega=omega;
    end
else
    omega1=(1+domega)*mu;
    omega2=(1-domega)*mu;
    omega = [omega1; omega2];
    if NA
        Omega = 2*pi/ 0.05 / tfin; % modulation of K
        eps=min(K0/2, max(1, K0/5));
        epsO=[1, 1]*abs(domega)*fm*10;
        OmegaO=2*pi ./ [0.0753, 0.0523]/tfin ; % modulation of omega
        %
        Oeff(1:Nt,1:2)=repmat(omega', Nt,1) + repmat(epsO,Nt,1).*sin(t'*OmegaO);
        Keff(1:Nt)=K_cnctm(1,2)+eps*sin(Omega*t);
    end
end
varstrKM.h=h;
varstrKM.Nds=Nds;
varstrKM.t=t;
varstrKM.mu=mu;
varstrKM.cnt=round(cnt0/dt);
varstrKM.N=N;
varstrKM.N2=N/2;
varstrKM.K_cnctm=K_cnctm;
varstrKM.idx=idx;
varstrKM.idx1=idx1;
%----
while i<Nh
    if N==2
        if NA
            varstrKM.omega=Oeff(ceil((i+1)/Nds),:)';
            K_cnctm(K_cnctm~=0)=Keff(ceil((i+1)/Nds));
            varstrKM.K_cnctm=K_cnctm;
        else
            varstrKM.omega=omega;
        end
    end
    eta=sqrt(2*h*D)*randn(N,1);
    chck=ceil((i+1)/Nds); 
    if min(taul(:))==0
        [theta, thds(:,chck+1), i]    = KMcnctmHt0PlosCB(theta, eta, i, chck, varstrKM);
    else
        [theta, thds(:,chck+1), i, ~] = KMcnctmHfPlosCB(theta,  eta, i, chck, varstrKM); 
    end
end
%% ---------------
% still to be able to print in case of breaking the simulation
[~,I]=find(isnan(thds), 1, 'first');
if isempty(I),
    I=length(thds(1,:))+1;
else
    if I<length(t)-1
        warning('NaN detected|')
        disp(['I=', int2str(I), ',  t=', num2str(t(I))])
        pause;
    end
end;
t(I:end)=[];
Nt=length(t);
tfin=floor(t(I-1));
thds(:,I:end,:)=[];
if NA
    Oeff=Oeff(1:Nt,:);
    Keff=Keff(1:Nt);
end
%% -------------
if N>2
    fname=[tit int2str(N) 'f' int2str(fm)  't'  int2str(tfin) 't' num2str(taul)...
        'K' num2str(K0) 'D' num2str(D) 'h' num2str(hpar) 'coh' num2str(incoh) 'v2']; 
else
    if NA
        fname=[int2str(N) 'f' int2str(fm) 'df' int2str(domega)  't'  int2str(tfin) 't' num2str(taul) 'K' num2str(K0) 'D' num2str(D) 'h' num2str(hpar) 'v2'];
    else
        fname=[int2str(N) 'f' int2str(fm) 'df' int2str(domega)  't'  int2str(tfin) 't' num2str(taul) 'K' num2str(K0) 'D' num2str(D) 'h' num2str(hpar) 'NAv2'];
    end
end
fname=regexprep(fname, '\s+', '');
fname=regexprep(fname,{'\.','txt'},{'','.txt'});
if N==2
    WL=8*2*pi/mean(mu);
    overlap=0.8; % overlapping of windows
    if NA
        [R,xR,surr,surr2,thsync2,thsyncplv2,thsync,thsyncplv,Dthmean,Dthvar,Dplvmean,Dplvvar,Dthmean2,Dthvar2,Dplvmean2,Dplvvar2,dth, PLVs, phiplv]=...
            plvstatPlosCB(thds, WL, overlap, Oeff);
%         printDph_ch([fname, int2str(100*domega)], 1, 2, taul, taul, 3, squeeze(dth(1,2,:)), surr, surr2(1,2),  PLVs(1,2), squeeze(R(1,2,:)), squeeze(phiplv(1,2,:)), xR, ...
%             squeeze(thsync(1,2,:)), squeeze(thsync2(1,2,:)), squeeze(thsyncplv(1,2,:)), squeeze(thsyncplv2(1,2,:)),...
%             0, '2 osc',[],[], Oeff, Keff(:));
    else
        [R,xR,surr,surr2,thsync2,thsyncplv2,thsync,thsyncplv,Dthmean,Dthvar,Dplvmean,Dplvvar,Dthmean2,Dthvar2,Dplvmean2,Dplvvar2,dth, PLVs, phiplv]=...
            plvstat(thds, WL, overlap);
%         printDph_ch([fname, int2str(100*domega)], 1, 2, taul, taul, 1, squeeze(dth(1,2,:)), surr, surr2(1,2),  PLVs(1,2), squeeze(R(1,2,:)), squeeze(phiplv(1,2,:)), xR, ...
%             squeeze(thsync(1,2,:)), squeeze(thsync2(1,2,:)), squeeze(thsyncplv(1,2,:)), squeeze(thsyncplv2(1,2,:)),...
%             0, '2 osc');
    end
    return
    %%
else
    %%
    par.h=h; par.dt=dt; par.tfin=tfin; par.D=D; par.fm=fm; par.hpar=hpar; par.N=N; par.Keff=K; par.K=K0; par.Nfmax=Nfmax;
    papar.taul=taul; 
    %
    %% print time-series of order parameters
    %[zs, zs1, zs2]=printFTevol(thds, fname, psave, 0.01, rl, NaN);
    zs = mean(exp(1j*thds),1);
    zs1 = mean(exp(1j*thds(1:N/2,:)));
    zs2 = mean(exp(1j*thds(N/2+1:N,:)));
    Omega = diff(unwrap(angle(zs)))/dt/2/pi;
    Omega=mean(Omega(round(Nt/2):end));
    Omega1 = diff(unwrap(angle(zs1)))/dt/2/pi;
    zs1=zs1(:);
    zs2=zs2(:);
    zmean=mean([abs(zs1(floor(3*Nt/5):Nt-1)); abs(zs2(floor(3*Nt/5):Nt-1))]);
    phi=angle([zs1(:),zs2(:)]);
    phi4=phi(floor(4*Nt/5):Nt-1,:);
    phasediff4 = diff(unwrap(phi4))/dt; 
    meanphase4=mean(mean(phasediff4(floor(0.5*length(phasediff4)):end-1,:)));
    delt=meanphase4*(taul(2)-taul(1))/2; % \Omega \Delta \tau
    beta=meanphase4*mean(taul); % \tilde{\Omega} \tau
    disp(['O*tau2= ', num2str(meanphase4*taul(2)/pi/2, '%6.2f'), ' pi, O*Dtau= ', num2str(delt, '%6.2f'), ...
        ' pi , O*tau1= ', num2str(meanphase4*taul(1)/pi/2, '%6.2f'), ' pi'])
%     %% plot entrained oscillators
%     if freqnoise
%         phi4=phi(floor(2*Nt/5):Nt-1,:);
%     end
%     theta=NaN(Nt,N/2,2);
%     theta(:,:,1)=thds(1:N/2,:)';
%     theta(:,:,2)=thds(1+N/2:end,:)';
%     thetan=mod(rem(theta+pi, pi*2)+pi*2, pi*2)-pi;    
%     phasediff2 = diff(unwrap(phi))/dt;
%     phasediff4 = diff(unwrap(phi4))/dt;
%     %%
%     if ~freqnoise        
%         figname=['randK'*Klog, tit, int2str(N), 'p' int2str(round(10*p1)) 'tau' num2str(taul(1)) '_' num2str(taul(2)) ...
%             'K' num2str(K0) 'O' int2str(round(mu(1)/pi/2)) 'g' int2str(D) ];
%         figname = regexprep(figname,'[^a-zA-Z0-9]','');
%         printOi_ph(thds, Arnd, zs1(:), zs2(:), mu, gamma, omega1, omega2, taul, dt, t, K0, ki, figname);        
%     else
%         relativeph=1; % if 1 then the phases in the right hemisphere are plot relative to those in the left
%         printsep=2; % print figures separately
%        printKi_ph(thds, zs1(:), zs2(:), mu, D, ki(:)', taul, dt, t, K0, fname, relativeph, printsep); % for random tau probably gives wrong results.
%     end
end

