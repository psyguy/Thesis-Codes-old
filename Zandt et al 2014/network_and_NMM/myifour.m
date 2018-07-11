function [y,t]=myifour(Y,Df,N)
% de sample frequentie  (dus niet correct)
if nargin < 3; N = length(Y); end
if length(Df) > 1; Df = Df(2)-Df(1); end
Dt = 1/(Df*N);
t = (0:N-1)*Dt;

Y = Y/Dt; % to get the right amplitude
y = ifft(ifftshift(Y),N);




