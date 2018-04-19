function LabelsF = updateLabels (Labels1, Labels2)

% UPDATELABELS: Update labels
%
% USAGE: LabelsF = updateLabels (Labels1, Labels2)
%
% Example:
%
% Labels1: 1 1 1   2 2  3  4 4  5 5 5
%
% Labels2:    1     1   2   1     2
%
% LabelsF: 1 1 1   1 1  2  1 1  2 2 2
%
%
%
% Author: Guiomar Niso, 2014


uLabels1 = unique(Labels1)'; 
% nLabels1 = length(uLabels1); % en teoria: nLabels1 == length(Labels2)

uLabels2 = unique(Labels2)';
nLabels2 = length(uLabels2);

LabelsF = zeros (size(Labels1));

for i=1:nLabels2

    iV2 = find( Labels2 == uLabels2(i) );

    for j=1:length(iV2)

        iV1 =  Labels1 == uLabels1(iV2(j));
        LabelsF(iV1) = uLabels2(i);

    end
end

end
