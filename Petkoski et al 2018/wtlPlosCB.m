%-------------------------------------------------------------------------%
%-                           Spase Petkoski                              -%
%- Theoretical Neuroscience Group -- Institute of Systemes Neuroscience  -%
%-                        2018, AMU Marseille                            -%
%-------------------------------------------------------------------------%
% The function wtlPlosCB returns weighted lengths of thelinks, where 
% weigths are integers. Links with weight w are counted w times. This
% allows calculating weighted means of the links within and between
% hemispheres. 
% It can also find indices of weights smaller than a threshold th and make 
% them to be NaN in the matrix of the delays (thus weaker links are 
% discarded)

function[tlw] = wtlPlosCB(tl,w,th)
  
    if th<0, error('th should be non-negative'); end
    if size(tl)~=size(w), error('matrices should have same size'); end
    tmp=tl; 
    tmp(w<=th)=NaN; 
    %----
    tlw=[];
    for i = find(~isnan(tmp))'
        tlw=[tlw, repmat(tmp(i),1,w(i))];   
    end
end
