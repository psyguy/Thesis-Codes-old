function output = NMMplotter(pMu,N)

index = 1; % to get only y
yinit = rand(1,4);%[0,0,0,0];
Y(1,:) = yinit;
for (i=1:N)
    p = normrnd(pMu,10);
%    p = pMu;
    Y(i+1,:) = grifau2006simple(p, Y(i,:));
end

output = Y(:,index);
end