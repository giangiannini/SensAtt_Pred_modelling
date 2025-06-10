function KL = KL_DIV(param1,param2)
% Alpha and Beta are vectors containing the alpha and beta parameters for
% each beta distribution that you want to calculate the difference of. 

%example
% param1 = [1 1]; 
% param2 = [1 1]; 

KL = gammaln(sum(param1)) - gammaln(sum(param2)) - sum(gammaln(param1)) + ...
sum(gammaln(param2)) + (param1 - param2) * (psi(param1) - psi(sum(param1)))';

end


