function p = BS_put(S, X, r, T, sigma, q)
  arguments
    S; X; r; T; sigma
    q=0
  end
  d1 = (log(S/X) + (r - q + sigma^2 / 2) * T) / (sigma * sqrt(T));
  d2 = d1 - sigma * sqrt(T);
  p = X * exp(-r * T) * normcdf(-d2) - ...
      S * exp(-q * T) * normcdf(-d1);
end