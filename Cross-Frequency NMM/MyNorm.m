function [ yn ] = MyNorm( y )
    yn = (y-mean(y))/std(y);
end

