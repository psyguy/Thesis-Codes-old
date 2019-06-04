function Set_ModelParam(hws,P,val)
%% load_system(MS.SimulinkFile)
    L = hws.data;
    [m] = arrayfun(@(x)find(strcmp(x.Name,P)),L,'uniformoutput',false);
    I = find(~cellfun(@isempty,m));
    L(I).Value.Value = val;

    
