function [ NRMSE ] = nrmse(data, data_ref)
% This function calculates NRMSE of data according to its reference
% 
% INPUTS:
%       data (N_by_M double): estimated data
%       data_ref (N_by_M double): reference/observed data
%
% OUTPUTS:
%       nrmse (1_by_1 double): nrmse

Xref = abs(data_ref);
X = reshape(abs(data), size(data_ref));
se = (Xref - X).^2;
mse = nanmean( se(:) );
rmse = sqrt( mse );
Xmean = nanmean(Xref(:));
NRMSE = rmse / Xmean;

end