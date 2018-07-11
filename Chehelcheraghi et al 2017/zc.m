function [n] = zc(Sig,th)
    nn2p=0;
    np2n=0;
    L = length(Sig);
    for i=1:L-1
       if (Sig(i)>th && Sig(i+1)<th) 
           np2n=np2n+1;
       end
       
       if (Sig(i)<th && Sig(i+1)>th) 
           nn2p=nn2p+1;
       end
        
    end

    n=nn2p+np2n;
end