function [TPR, FPR, thres] = check_roc (GT, RES)

% ROC curve
% 
% [TPR, FPR, thres] = check_roc (GT, RES)
%
% INPUT:
% - GT: Ground Truth
% - RES: Results to check
%
% OUTPUT:
% - TPR: True Positive Rate
% - FDR: False Positive Rate
% - thres: Number of points for the ROC curve
%
% Author: Guiomar Niso, 2014

% Number of points for the ROC curve
% thres = linspace( min(RES), max(RES), Nt);
thres = min(RES):max(RES);

Nt = length(thres);

TPR = zeros(1,Nt);
FPR = zeros(1,Nt);

for i = 1:Nt
    
    R = RES>=thres(Nt-i+1);
    
    TPR(i) = sum( R(GT==1) ) / sum( GT==1 );
    FPR(i) = sum( R(GT==0) ) / sum( GT==0 );
    
end

end
