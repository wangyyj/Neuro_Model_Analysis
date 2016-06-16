%Calc_R2
function [rsq coef pval] = calc_R2(x,y)
% http://www.mathworks.com/help/matlab/data_analysis/linear-regression.html
%Calculates R-squared value

if length(y)~=length(x)
    warning('Unequal numbers of X and Y')
end
%p(1) is the slope and p(2) is the intercept of the linear predictor
[coef stats] = robustfit(x,y); 

yfit =  coef(2) * x + coef(1);
% yresid = y - yfit;
% SSresid = sum(yresid.^2);
% SStotal = (length(y)-1) * var(y);
% rsq = 1 - SSresid/SStotal;
rsq = rsquare(y,yfit);

pval=stats.p(2);
