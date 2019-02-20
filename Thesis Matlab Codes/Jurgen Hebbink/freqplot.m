function [f,power]=freqplot(u,Fs)
%gebaseerd op voorbeeld help fft. Inputsignaal u, samplefrequentie Fs.
L=length(u);
NFFT=2^nextpow2(L);%punten voor transformatie
Y=fft(u,NFFT)/L;%fouriergetransformeerde
f=Fs/2*linspace(0,1,NFFT/2+1);%frequenties fouriertransformatie
power=2*abs(Y(1:NFFT/2+1));
%plot(f,power);

end