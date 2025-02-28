%{
function c = BS_call(S,X,r,T,sigma,q)
  arguments % default arg
    S; X; r; T; sigma
    q=0
  end
  d1 = (log(S/X) + (r-q+sigma^2/2)*T)/sigma/sqrt(T);
  d2 = d1 - sigma * sqrt(T);
  c = S.*(exp(-q*T)*normcdf(d1)) - X*exp(-r*T)*normcdf(d2);
end
%}

function cdo = BS_doCall(S0, X, r, T, sigma, q, H)
    lambda = (r - q + (sigma^2)/2) / sigma^2; 
    y = log(H.^2 ./ (S0.*X))./(sigma*(sqrt(T))) + lambda*sigma*sqrt(T);
    cdi = S0.*exp(-q*T).*(H ./ S0).^(2*lambda) .* normcdf(y) - X*exp(-r*T).*(H./S0).^(2*lambda - 2).*normcdf(y - sigma*(sqrt(T)));
    c = BS_call(S0, X, r, T, sigma, q);
    cdo = c - cdi;
end