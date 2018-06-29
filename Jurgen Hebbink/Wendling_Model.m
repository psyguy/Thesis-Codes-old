tic; %To measure computation time
% 
%% Simulation settings
Npart=2; % Number of parts for fourier transform including intial part
Tend=10*Npart;%Totalsimulationtime
NTST=Tend*10000;%NumberofTimeSteps
Nout=5;%Outputevery5thtimestep

deltat=Tend/NTST;%Timestepsize

%%Parametersettings
A=5;%Averageexcitorysynapticgain
B=25;%Averageslowinhibitorysynapticgain
G=15;%Averagefastinhibitorysynapticgain

a=100;%Inverseofexcitorytimeconstant
b=50;%Inverseofslowinhibitorytimeconstant
g=500;%Inverseoffastinhibitorytimeconstant

C=135;%Connectivityparameter
C1=C;
C2=0.8*C;
C3=0.25*C;
C4=0.25*C;
C5=0.3*C;
C6=0.1*C;
C7=0.8*C;

pf=90;%Meaninput
sd=2;%Standarddeviationinput

%%Initialcondititons
x=zeros(4,1);
y=zeros(4,1);

%%Initialisevariablesforcomputation
Im=eye(4,4);

Q1=-deltat*diag([a^2,a^2,b^2,g^2]);
Q2=Im-2*deltat*diag([a,a,b,g]);
Q3=deltat*diag([A*a,A*a,B*b,G*g]);


fd=[0;A*a*pf/C2;0;0]*deltat;
fs=[0;A*a/C2;0;0]*sqrt(deltat)*sd*randn(1,NTST);

P=[[0,C2,-C4,-C7];
[C1,0,0,0];
[C3,0,0,0];
[C5,0,-C6,0]];

u=P*x;
uout=zeros(4,ceil(NTST/Nout)+1);
uout(:,1)=u;

xn=zeros(4,1);
yn=zeros(4,1);

%%
t1=toc;
disp(['Initialisationtime:',num2str(t1)]);
disp('Startsimulation');

tel=0;
indexplot=2;
for i=1:NTST
xn=x+deltat*y;
yn=Q1*x+Q2*y+Q3*sigmf(u,[1,0])+fd+fs(:,i);

x=xn;
y=yn;
u=P*x;

%SavesuforoutputeveryNoutiterations
tel=tel+1;
    if tel==Nout
        tel=0;
        uout(:,indexplot)=u;
        indexplot=indexplot+1;
    end
end

t2=toc;
disp(['Simulationcomplete,duration:',num2str(t2-t1)]);

%%Postprocessing

%Plottimeseries
figure;
plot(0:deltat*Nout:Tend,uout(1,:));
xlim([0 Tend])
set(gca,'YDir','reverse');

%Calculatefrequencyspectrum
for j=1:Npart-1;
    [f,pow]=freqplot(uout(1,floor(NTST*j/(Npart*Nout)):floor(NTST*(j+1)/(Npart*ut)))-mean(uout(1,floor(NTST*j/(Npart*Nout)):floor(NTST*(j+1)/(Npart*ut)))),1/(deltat*Nout));
    if j>1;
        powtot=powtot+pow;
    else
        powtot=pow;
    end
end

%Displayfrequencyspectrum
figure;
plot(f,powtot/(Npart-1));
xlim([050])

toc; %Displaystotalcomputationtime