% calculate a discrete fourier tranform and use the sample time to create
% a frequency axis. Also, correct the amplitude for the sample time.
% [Y,f] = myfour(y,Dt,N); if y is 2D, then the fourier transform of all
% !columns! are calculated.

function [Y,f]=myfour(y,Dt,N)
if nargin < 3; N = length(y); end
    
if length(Dt)> 1; Dt = Dt(2)-Dt(1); end; % if complete time axis is given, assume evenly spaced samples

[a,b] = size(y);
if a == 1 || b == 1
    Y = fftshift(fft(y,N));
else
Y = fftshift(fft(y,N),1);
end
Y = Y*Dt; % to get the right amplitude

Df = 1/(Dt*N);
f = ((0:N-1)-N/2)*Df;