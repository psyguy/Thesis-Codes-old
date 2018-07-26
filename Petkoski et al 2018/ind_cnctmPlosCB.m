function[idx, ind, inddif, idxind1, idxind, idx1]=ind_cnctmPlosCB(theta, tau_cnctm)
% finding the indices of delayed phases
global N
NN=size(theta);
tau_m=NN(2);
%%%
[ind, inddif, idxind, idxind1]=deal(cell(N,1));
idx=NaN(N,N); idx1=idx;
for m=1:N
    % used to be tau_cnctm(:,m), though it doesn't matter
    ind{m}    = find(tau_cnctm(m,:));                                   % non-zero links for each node
    inddif{m} = setdiff(1:N,ind{m});                                    % zero links for each node
    idx(m,:)  = sub2ind(NN, 1:N,    tau_m-tau_cnctm(m,1:N));            % indices for delayed values of each phase (even those with 0 delays and/or weights)
    idxind1{m}= sub2ind(NN, ind{m}, tau_m-tau_cnctm(m,ind{m}) + 1);    % indices for delayed values +1, for phases with non-zero delays
    idx1(m,ind{m}) = idxind1{m}; % indices for delayed values +1, for phases with non-zero delays; the others are just delayed
end
idx1(isnan(idx1))=idx(isnan(idx1));



