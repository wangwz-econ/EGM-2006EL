function [EUP] = fun_GothVP(a)
% GothVP is defined as Gothic V prime in the paper 
%        (end of period marginal value). 

global beta rho alpha dep g R W u_D_inv u_D permval permprob

% EUP is the sum of marginal utilities UP weighted by probabilities

EUP = zeros(size(a)); 

for j = 1: length(permval)

    CapPsi = permval(j);
    
    k = a .* dep ./ (g .* CapPsi);
    EUP = EUP + dep * beta * R(k) .* ... 
        u_D( g .* CapPsi .* fun_Cnextp(R(k).*k + W(k)) ) .* permprob(j);

end

izero = (a==0.0);
EUP(izero) = Inf;
