%-------------------------------------------------------------------------%
%-                           Spase Petkoski                              -%
%- Theoretical Neuroscience Group -- Institute of Systemes Neuroscience  -%
%-                        2018, AMU Marseille                            -%
%-------------------------------------------------------------------------%
function[R,xR,surr1,surr2,thsync2,thsyncplv2,thsync,thsyncplv,Dthmean,Dthvar,...
    Dplvmean,Dplvvar,Dthmean2,Dthvar2,Dplvmean2,Dplvvar2, dth, PLVs, phiplv]=...
    plvstatPlosCB(thds, WL, overlap, varargin)
%---
% function that first generates the surrogates according to the two
% procedures: (i) uncoupled oscillators with the same parameters as in the
% time series, and (ii) shuffled time series; and then calculates the
% levels of significance for the coherence, as described in the manuscript,
% and then finally calculates cPLV and henceforth |PLV| and angle(PLV), as
% well as other metrics derived from the phase relationships.
% The length of the time-windows and the overlap is the same for the
% original and for the surrogate time-serie, and for each time-window a 
%
% thds - matrix of phases from N regions at Nt time points
% WL - window length (usualy 8 to 10 periods)
% overlap -  determines how much two consecutive windows overlap; 
%            overlap =  (0, 1];  (usualy 0.75); 1 for no overlap.
% omegaN - natural frequencies for each oscillator (can be time-varying)
% R - PLV (absolute value) for each pair at each time-window
% PLVs - PLV calculated for the whole-time series
% xR - middle points of each window
% surr1 - first level of statistical significnce (uncoupled oscillators)
% surr2 - second level of statistical significnce (shuffled phases)
% dth - pair-wise difference of angles
% phiplv - angle of PLV in each time-window in (-pi, pi)
% thsync - has the values of dth when R>surr1, and it is NaN otherwise
% thsync2 - has the values of dth when R>surr2, and it is NaN otherwise
% thsyncplv - has the values of phiplv when R>surr1, and NaN otherwise
% thsyncplv2 - has the values of phiplv when R>surr2, and NaN otherwise
% Dthmean - circular mean of thsync 
% Dthvar - circular variance of thsync 
% Dthmean2 - circular mean of thsync2
% Dthvar2 - circular variance of thsync2
% Dplvmean - circular mean of thsyncplv 
% Dplvvar - circular variance of thsyncplv 
% Dplvmean2 - circular mean of thsyncplv2 
% Dplvvar2 - circular variance of thsyncplv2 

global dt D omega mu
%
[N, Nt]=size(thds);
pi2=2*pi;
% omega added for 2 oscillators when there is NA
if ~isempty(omega)
    omegaN=repmat(omega(:)',Nt,1);
else
    omegaN=mu*ones(Nt,2);
end
if nargin==4 && ~isempty(varargin{1}),
    omegaN=varargin{1}; % Oeff
end
% calculate PLV only after the first 10 seconds to avoid transient effect.
T=round(10/dt); % transition points, for t<10s.
Nsurr=100; % number of surrogates time-series
wnd = round(WL/dt); % number of data points in one window
ds  = round(wnd*overlap);
[R, phiplv]=deal(NaN(      N,N,floor((Nt-wnd)/ds)+1 ));
surr2=NaN(N,N);
%%
for m=1:N
    for n=m+1:N
        [R(m,n,:),  xR, phiplv(m,n,:)] = plvPlosCB(thds(m,:), thds(n,:), dt, WL, overlap); % theta_n-theta_m
        R(n,m,:)=R(m,n,:);
    end
end
%%
% surrogates; noise driven, uncoupled oscillators
% for N=2, frequencies can differ
DD=sqrt(2*D*dt);
ccnt=0;
Rp1=  NaN(Nsurr,    floor((Nt-wnd)/ds)+1 ); % although for Rp1 all are the same
for ss=1:Nsurr
    thnc=NaN(2,Nt);
    thnc(:,1)=rand(2,1)*2*pi;
    eta=randn(2,Nt-1);
    for m=1:Nt-1
        thnc( :,m+1)=thnc(:,m) + DD*eta(:,m) + omegaN(m,:)'*dt; % Euler scheme; seems to be the same as Heun for uncoupled case!
    end
    [Rp1(ss,:), ~, ~] =  plvPlosCB(thnc(1,:), thnc(2,:), dt, WL, overlap);
end
surr1=prctile(prctile(Rp1', 95), 95);

% surrogates, shuffled time series
for m=1:N
    for n=m+1:N
        Rp2=NaN(Nsurr,    floor((Nt-wnd)/ds)+1 ); % although for Rp1 all are the same
        for ss=1:Nsurr
            [Rp2(ss,:), ~, ~] = plvPlosCB(thds(m,:), thds(n,randperm(Nt)), dt, WL, overlap);
        end
        surr2(m,n)=prctile(squeeze(max(Rp2,[], 2)),95);
        surr2(n,m)=surr2(m,n);
        ccnt=ccnt+1;
        if ~mod(ccnt,400),
            disp(['plv, ', num2str(ccnt), '; out of ', num2str(N*(N-1)/2), '; ',num2str(100*ccnt/(N*(N-1)/2), '%3.1f'), '%'])
        end
    end
end
clear Rp*
%% --
dth=NaN(N,N,Nt);
for m=1:N
    for n=m+1:N
        dth(m,n,:)=mod(rem(thds(n,:)-thds(m,:)+pi, pi2)+pi2, pi2)-pi; % delta theta dth(n,m,:)=-dth(m,n,:);
    end
end
thsync2=NaN(N,N,Nt); thsync=thsync2;
thsyncplv=NaN(N,N,floor((Nt-wnd)/ds)+1); thsyncplv2=thsyncplv;
TMP=NaN(N,N);
TMP(logical(eye(size(TMP))))=0; % temporary matrix with zero diagonal.
[Dthmean,Dthvar,Dplvmean,Dplvvar,Dthmean2,Dthvar2,Dplvmean2,Dplvvar2,PLVs]=deal(TMP);
clear TMP;
ccnt=0;
for m=1:N
    for n=m+1:N
        for ii = 0:floor((Nt-wnd)/ds)
            ids=ii*ds;
            if R(m,n,ii+1) > surr1 && ids>T
                L=ids+1:ids+wnd;
                thsync(m,n,L)=dth(m,n,L);
                thsync(n,m,L)=thsync(m,n,L);
                thsyncplv(m,n,ii+1)=phiplv(m,n,ii+1);
                thsyncplv(n,m,ii+1)=thsyncplv(m,n,ii+1);
            end
            if R(m,n,ii+1) > surr2(m,n) && ids>T
                L=ids+1:ids+wnd;
                thsync2(m,n,L)=dth(m,n,L);
                thsync2(n,m,L)=thsync2(m,n,L);
                thsyncplv2(m,n,ii+1)=phiplv(m,n,ii+1);
                thsyncplv2(n,m,ii+1)=thsyncplv2(m,n,ii+1);
            end
        end
        PLVs(m,n)=abs(KOP(squeeze(dth(m,n,T:end)))); PLVs(n,m)=PLVs(m,n);
        % mean and var of the phase lags, stationarity assumed.
        thtmp=squeeze(thsync(m,n,:));
        thtmp=thtmp(~isnan(thtmp));
        if ~isempty(thtmp)
            Ztmp=KOP(thtmp);
            Dthmean(m,n)=angle(Ztmp);
            Dthmean(n,m)=-Dthmean(m,n);
            Dthvar(m,n) =1-abs(Ztmp);
            Dthvar(n,m) =Dthvar(m,n);
        end
        %
        thtmp=squeeze(thsyncplv(m,n,:));
        thtmp=thtmp(~isnan(thtmp));
        if ~isempty(thtmp)
            Ztmp=KOP(thtmp);
            Dplvmean(m,n)=angle(Ztmp);
            Dplvmean(n,m)=-Dplvmean(m,n);
            Dplvvar(m,n) =1-abs(Ztmp);
            Dplvvar(n,m)=Dplvvar(m,n);
        end
        %
        thtmp=squeeze(thsync2(m,n,:));
        thtmp=thtmp(~isnan(thtmp));
        if ~isempty(thtmp)
            Ztmp=KOP(thtmp);
            Dthmean2(m,n)=angle(Ztmp);
            Dthmean2(n,m)=-Dthmean2(m,n);
            Dthvar2(m,n) =1-abs(Ztmp);
            Dthvar2(n,m)=Dthvar2(m,n);
        end
        %
        thtmp=squeeze(thsyncplv2(m,n,:));
        thtmp=thtmp(~isnan(thtmp));
        if ~isempty(thtmp)
            Ztmp=KOP(thtmp);
            Dplvmean2(m,n)=angle(Ztmp);
            Dplvmean2(n,m)=-Dplvmean2(m,n);
            Dplvvar2(m,n) =1-abs(Ztmp);
            Dplvvar2(n,m) =Dplvvar2(m,n);
        end
        
        ccnt=ccnt+1;
        if ~mod(ccnt,500),
            disp(['significance, ', num2str(ccnt), '; out of ', num2str(N*(N-1)/2), '; ',num2str(100*ccnt/(N*(N-1)/2), '%3.1f'), '%'])
        end
    end
end

end
