close all;
fl = 15;
fh = 0;
fs = 10000;
LL = 1.5;
offset = 2; % seconds
t = (0:1/fs:LL);


%phase_regions(so_1,out_f_2,1,0,fs,0);
figure('units','normalized','outerposition',[0 0 1 1]); %%% Figure 3
%suptitle('Slow and Output, Filtered');
subplot(1,2,1);
title('Node 1','FontSize',20);
hold on;
ylim([-2 8]);
xlim([0 2]);
phase_regions(Node_1,Node_1,Node_2,Node_2,500,0,offset,LL,fs,fl,fh,4);
legend({'Output Lower-Frequency Band','Output Upper-Frequency Band'},'FontSize',12,'FontWeight','bold');
legend boxoff

subplot(1,2,2);
title('Node 2','FontSize',20);
hold on;
ylim([-2 8]);
phase_regions(Node_2,Node_2,Node_1,Node_1,500,0,offset,LL,fs,fl,fh,4);
legend({'Output Lower-Frequency Band','Output Upper-Frequency Band'},'FontSize',12,'FontWeight','bold');
legend boxoff
%set(lgnd,'color','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

