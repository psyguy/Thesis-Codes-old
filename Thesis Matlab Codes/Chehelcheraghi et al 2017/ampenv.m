function [env] = ampenv(y,g)

% Reads in the sound file, into a big array called y.


% Normalize y; that is, scale all values to its maximum. Note how simple it is 
% to do this in MatLab.

%y = y/max(abs(y));

% If you want to hear the sound after you read it in, uncomment this next line. 
% Sound just plays a file to the speaker.

%sound(y, 44100);

% Initialize peak and average arrays to be the same as the signal. We�re 
% keeping three big arrays, all of the same length.

averageenv = y;						
maxenv = y;

% Set the window lengths directly, in number of samples. Generally, the average 
% should be smaller than the peak, since it will tend to flatten out if it gets 
% too big. These two numbers can be played with to vary what the pictures will 
% look like. Note that we keep two different windowsizes so that we can have 
% different kinds of resolution for peak and average charts.

averagewindowsize = 16;						
peakwindowsize = 24;

% Go through the input signal, taking an average of the previous, 
% averagewindowsize number of samples, and store that in the current place in
% the average array. We do this by having a loop that starts at the end of the
% first window and goes to the end of the sound file, indicated by 
% length(averageenv). The MatLab command sum takes some range of a 
% vector, and gives you back the sum.

for k = averagewindowsize:length(averageenv)							
runningsum = sum(abs(y(k-averagewindowsize+1:k)));							
averageenv(k) = runningsum /averagewindowsize ;							
end

% Go through the input signal, taking the maximum value of the previous 
% peakwindowsize number of samples. We do this in the same way as we did the 
% average, but now we use max instead of sum.

for k = peakwindowsize:length(maxenv)							
maxenv(k) = max(y(k-peakwindowsize+1:k));							
end

% Plot the three "signals": the original in the color cyan, the 
% peak in magenta, and the average in blue. The command hold on invokes the 
% ghost of Sam and Dave to allow us to overprint on the current graph, so that 
% we can show all three lines together. The rest of this is just MatLab 
% graphics instructions (titles, labels, axes, etc.). The colors are the 
% letters in single quotes.

% figure;
% plot(y, 'c')						
% hold on
% plot (maxenv, 'm')						
% hold on
% plot(averageenv, 'b')
% title('Monochord: Signal, Average Signal Envelope, Peak Signal Envelope')
% xlabel('sample number')
% ylabel('amplitude (-1 to 1)')
% legend('original signal', 'peak envelope', 'running average envelope', 0)

if (g==1)
    env=maxenv;
else
    env=averageenv;
end