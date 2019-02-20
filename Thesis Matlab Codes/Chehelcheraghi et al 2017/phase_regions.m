function [n,nf] = phase_regions(SO1,FO1,SO2,FO2,delay,th,offset,LL,fs,fl,fh,flag)

%     SO1 = MyNorm(SO1);
%     FO1 = MyNorm(FO1);
    FO1 = FO1(delay:end-1);
    SO1 = SO1(1:end-delay);
    m = 0;
   
    if (flag==1)
        [FO1] = MyNorm(eegfilt(FO1',fs,fl,fh));
        [FO2] = MyNorm(eegfilt(FO2',fs,fl,fh));
        [SO1] = MyNorm(eegfilt(SO1',fs,fh,fl));
        [SO2] = MyNorm(eegfilt(SO2',fs,fh,fl));
    elseif (flag == 2)
        %[b,a]= cheby1(8, 0.5, 0.01, 'high');
        [b,a] =butter(4, 0.0005, 'high') ;    
        [FO1] = filter(b,a,FO1);  
    elseif (flag == 3)
        FO1 = FO1-SO1;
    elseif (flag == 4)
        [FO1] = eegfilt(FO1',fs,fl,fh);
        [FO2] = eegfilt(FO2',fs,fl,fh);
        [SO1] = eegfilt(SO1',fs,fh,fl);
        [SO2] = eegfilt(SO2',fs,fh,fl);
    end
    
%     FO1 = MyNorm(FO1(offset*fs:(offset+LL)*fs));
%     SO1 = MyNorm(SO1(offset*fs:(offset+LL)*fs));

    [FO1, ~, norm_f] = MyMutualNorm(FO1(offset*fs:(offset+LL)*fs),FO2(offset*fs:(offset+LL)*fs));
    [SO1, ~, norm_s] = MyMutualNorm(SO1(offset*fs:(offset+LL)*fs),SO2(offset*fs:(offset+LL)*fs));
    %Env1 = ampenv(FO1,1);
    nn2p=0;
    np2n=0;
    L = length(SO1);
    t = linspace(0,L/fs,L);
    %figure;
    
    hold on;
    plot(t,SO1,'LineWidth',1.5);
    plot(t,FO1+3,'color','k','LineWidth',1.5);
    %plot(t,Env1+4.5,'color','r','LineWidth',1.5);
    
    ss1 = 1;
    ee1 = 1;
    ll = -1.5;
    hh = 4.5;
    
    for i=1:L-1
       if (SO1(i)>th && SO1(i+1)<th) 
           %% Detect Top to Down
           ss1 = i;
           np2n=np2n+1;
           %vline(i/fs,'r',num2str(floor((nf*fs)/(2*abs(ss1-ee1-2*m)))));
           patch([ss1/fs ss1/fs ee1/fs ee1/fs],[ll hh hh ll],'yellow','FaceAlpha',.1,'EdgeColor','none');
           nf = zc(FO1(ee1+m:ss1-m),th); 
           if (ee1>1)
            text((ss1+ee1-300)/(2*fs),2,num2str(floor((nf*fs)/(2*abs(ss1-ee1-2*m)))),'FontSize',16,'FontWeight','bold');
           end
       end
       
       if (SO1(i)<th && SO1(i+1)>th)
           %% Detect Down to Top
           ee1 = i;
           nn2p=nn2p+1;
           %vline(i/fs,'k',num2str(floor((nf*fs)/(2*abs(ss1-ee1-2*m)))));
           patch([ss1/fs ss1/fs ee1/fs ee1/fs],[ll hh hh ll],'green','FaceAlpha',.1,'EdgeColor','none');
           nf = zc(FO1(ss1+m:ee1-m),th); 
           if (ss1>1)
            text((ss1+ee1-300)/(2*fs),2,num2str(floor((nf*fs)/(2*abs(ss1-ee1-2*m)))),'FontSize',16,'FontWeight','bold');
           end
       end
       

       
    end
    
    if (SO1(L)>th)
        patch([ee1/fs ee1/fs L/fs L/fs],[ll hh hh ll],'yellow','FaceAlpha',.1,'EdgeColor','none');
    end
    
    if (SO1(L)<th)
        patch([ss1/fs ss1/fs L/fs L/fs],[ll hh hh ll],'green','FaceAlpha',.1,'EdgeColor','none');
    end
    n=nn2p+np2n;
    xlim([0 LL]);
    set(gca,'YTick',[-1.5 -1 0 1 2 3 4 4.5]);
    set(gca,'YTickLabel',{'',num2str(-norm_s),'0',num2str(norm_s),num2str(-norm_f),'0',num2str(norm_f),''});
    ylim([ll hh]);
    hline(1.5,'k');
    ylabel('Voltage [$mV$]','fontsize',20,'interpreter','latex','FontName', 'Times New Roman-Normal');
    xlabel('time [$sec$]','fontsize',20,'interpreter','latex','FontName', 'Times New Roman-Normal');
end