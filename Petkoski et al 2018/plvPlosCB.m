%-------------------------------------------------------------------------%
%-                           Spase Petkoski                              -%
%- Theoretical Neuroscience Group -- Institute of Systemes Neuroscience  -%
%-                        2018, AMU Marseille                            -%
%-------------------------------------------------------------------------%
%  Function for calcuulating phase locking values.
%       
%      Usage: [R, x, phi] = plv(th1,th2,dt,WL,k);
%
%  th1, th2 - phases
% dt - sampling interval --- 1/fs
% WL - window length (usualy 8 to 10 periods)
% k -  determines how much two consecutive windows overlap; k =  (0, 1]; 
% (usualy 0.75); 1 for no overlap.
% R - PLV (absolute value)
% phi - angle of PLV in (-pi, pi)
% x - middle points of each window
function [R, x, phi] = plvPlosCB(th1,th2,dt,WL,k)
if length(th1)~=length(th2), error('phases should have the same length'); end;
wnd = round(WL/dt); % number of data points in one window
ds  = round(wnd*k);
Z = zeros(floor((length(th1)-wnd)/ds)+1,1);

th1=th1(:);
th2=th2(:);
for i = 0:floor((length(th1)-wnd)/ds)
    ids=i*ds;
    Ind=ids+1 : ids+wnd;
    Z(i+1)=mean(exp(1j* ( th2(Ind)-th1(Ind) ) ));    
end
x = wnd/2 : ds : length(th1)-wnd/2;
x = x*dt;
R=abs(Z);
phi=angle(Z);
