function [ yn1, yn2, norm ] = MyMutualNorm(y1,y2)
    
%     yn1 = (y1-mean(y1))/sqrt(std(y1)*std(y2));
%     yn2 = (y2-mean(y2))/sqrt(std(y1)*std(y2));
%     norm = sqrt(std(y1)*std(y2));

     norm = max(max(y1-mean(y1)),max(y2-mean(y2)));
     yn1 = (y1-mean(y1))/norm;
     yn2 = (y2-mean(y2))/norm;
     
end