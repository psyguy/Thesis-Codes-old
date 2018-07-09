tic; %To measure computation time
%% Simulation settings
Npart=3;
Tend=3*Npart;                %Simulation time
NTST=Tend*10000;        %Number of Time Steps
Nout=5;                 %Output every 5th timestep

deltat=Tend/NTST;       %Timestep

%% Parameter settings
A=3.25;                 %Average excitory synaptic gain
B=22;                   %Average slow inhibitory synaptic gain
G=10;                   %Average fast inhibitory synaptic gain
%A=7, B=30, G=150 or A=5, B=30, G=70 for pSWD case


a=100;                  %Inverse of excitory time constant
b=50;                   %Inverse of slow inhibitory time constant
g=500;                  %Inverse of fast inhibitory time constant

C=135;                  %Connectivity parameter
C1=C;
C2=0.8*C;
C3=0.25*C;
C4=0.25*C;
C5=0.3*C;
C6=0.1*C;
C7=0.8*C;

pf=90;%90;                  %Mean input
sd=2;                  %Standard deviation input

%% Initial condititons
x=zeros(4,1);
y=zeros(4,1);

%% Initialise variables for computation
 Im=eye(4,4);
% Imdt=deltat*Im;

Q1=-deltat*diag([a^2,a^2,b^2,g^2]);
Q2=Im-2*deltat*diag([a,a,b,g]);
Q3=deltat*diag([A*a,A*a,B*b,G*g]);

fd=[0;A*a*pf/C2;0;0]*deltat;

fs=[0;A*a/C2;0;0]*sqrt(deltat)*sd*randn(1,NTST);

P=[[0, C2, -C4,-C7];
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
disp(['Initialisation time: ',num2str(t1)]);
disp('Start simulation');

tel=0;
indexplot=2;
for i=1:NTST
    xn=x+deltat*y;
    yn=Q1*x+Q2*y+Q3*sigm(u)+fd+fs(:,i);

    x=xn;
    y=yn;
    u=P*x;
    
    %Save u for output every Nout iterations
    tel=tel+1;
    if tel==Nout
        tel=0;
        uout(:,indexplot)=u;
        indexplot=indexplot+1;
    end
end

t2=toc;
disp(['Simulation complete, duration: ',num2str(t2-t1)]);
%% Post processing
for k=1:4
    figure;
    plot(0:deltat*Nout:Tend,uout(k,:));
    xlim([0 Tend])
end
set(gca,'YDir','reverse');

for j=1:Npart-1;
    j
    [f,pow]=freqplot(uout(1,floor(NTST*j/(Npart*Nout)):floor(NTST*(j+1)/(Npart*Nout)))-mean(uout(1,floor(NTST*j/(Npart*Nout)):floor(NTST*(j+1)/(Npart*Nout)))),1/(deltat*Nout));
    if j>1;
        powtot=powtot+pow;
    else
        powtot=pow;
    end
end
figure;
plot(f,powtot/(Npart-1));
xlim([0 50])
toc;