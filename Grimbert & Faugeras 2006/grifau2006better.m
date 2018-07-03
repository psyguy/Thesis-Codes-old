function Y = grifau2006better(varargin)
% p, n, Yprev

%% defining parameters

A=3.25;
B=22;
a=100;
b=50;
C=135;
c1=C;
c2=0.8*C;
c3=0.25*C;
c4=0.25*C;

aa = A/a;
bb = B/b;

% Parameters of the sigmoid
vmax = 5;
r = 0.56;
v0 = 6;
parsigm =[v0,vmax,r];

%% Init

defaults = {0,1,rand(1,4)};
defaults(1:nargin) = varargin;

p = cell2mat(defaults(1));
n = cell2mat(defaults(2));
Yprev = cell2mat(defaults(3));

yp = Yprev(1);
y0p = Yprev(2);
y1p = Yprev(3);
y2p = Yprev(4);

%% solve

for (i=1:n)
    
    y=yp; y0=y0p; y1=y1p; y2=y2p;
    
    y0n = aa*sigm(y1-y0,parsigm);
    y1n = aa*(p+c2*sigm(c1*y0,parsigm));
    y2n = bb*c4*sigm(c3*y0,parsigm);
    yn = aa*p + aa*c2*sigm(aa*c1*sigm(y,parsigm),parsigm) - bb*c4*sigm(aa*c3*sigm(y,parsigm),parsigm);
    
    Y(i,1)=yn; Y(i,2)=y0n; Y(i,3)=y1n; Y(i,4)=y2n;
    yp=yn; y0p=y0n; y1p=y1n; y2p=y2n;
    
end