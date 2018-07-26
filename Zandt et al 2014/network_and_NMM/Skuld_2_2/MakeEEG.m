function EEG = MakeEEG(Im, PosData)
% This fucntion determines the field generated by the SPyr cells

nCell = size(Im,2);

Depth = PosData(1:nCell,3);

Contribution = Im;

for iCell = 1:nCell
    Contribution(:,iCell) = Contribution(:,iCell)/Depth(iCell);
end

EEG = -sum(Contribution,2);