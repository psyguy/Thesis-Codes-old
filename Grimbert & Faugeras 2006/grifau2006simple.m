function Y = grifau2006simple(p, Yprev)

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
y = Yprev(1);
y0 = Yprev(2);
y1 = Yprev(3);
y2 = Yprev(4);

%% solve
y0 = aa*sigm(y1-y0,parsigm);
y1 = aa*(p+c2*sigm(c1*y0,parsigm));
y2 = bb*c4*sigm(c3*y0,parsigm);
y = aa*p + aa*c2*sigm(aa*c1*sigm(y,parsigm),parsigm) -bb*c4*sigm(aa*c3*sigm(y,parsigm),parsigm);

Y = [y, y0, y1, y2];
end