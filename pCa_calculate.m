function [hillCoeff ec50]=pCa_calculate(X_pCa,Y_pCa)
%   Inputs
%       Two vectors of the same length, the first containing the concentration and
%       the second the force. If doses of 0 are contained in the data,
%       these are used as control values and the data are normalised by
%       their mean value.

%deal with 0 dosage by using it to normalise the results.
normalised=0;
if (sum(X_pCa(:)==0)>0)
    %compute mean control response
    controlResponse=mean(Y_pCa(X_pCa==0));
    %remove controls from dose/response curve
    Y_pCa=Y_pCa(X_pCa~=0)/controlResponse;
    X_pCa=X_pCa(X_pCa~=0);
    normalised=1;
end

%hill equation sigmoid
sigmoid=@(beta,x)beta(2)+(beta(1)-beta(2))./(1+(x/beta(3)).^beta(4));
% cf_23(x) = D+(A-D)/(1+(x/C)^B)
% A is the lower asymptote so guess it with min(y)
% B is the Hill's slope so guess it with the slope of the line between first and last point.
% C is the inflection point (the concentration of analyte where you have
% half of the max response) so guess it finding the concentration whose
% response is nearest to the mid response.
% D is the upper asymptote so guess it with max(y)

%calculate some rough guesses for initial parameters
minResponse=min(Y_pCa);
maxResponse=max(Y_pCa);
slope=(Y_pCa(end)-Y_pCa(1))/(X_pCa(end)-X_pCa(1));
[~,Idx]=min(abs((Y_pCa-((max(Y_pCa)-min(Y_pCa))/2))));

%fit the curve and compute the values
% [coeffs,r,J]=nlinfit(dose,response,sigmoid,[minResponse maxResponse midResponse 1]);
[coeffs,r,J]=nlinfit(X_pCa,Y_pCa,sigmoid,[minResponse maxResponse X_pCa(Idx) sign(slope)]);
ec50=coeffs(3);
hillCoeff=coeffs(4);

