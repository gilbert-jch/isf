function [V] = data_coxeter(n,fout,verb)

%%
% [V] = data_coxeter(n,fout,verb)
%
% Provides the Coxeter arrangement in Rn defined by the hyperplanes
% H(i,j) := {x in Rn: x_i = x_j}. It is known that there is n! chambers.

  fprintf('\nCoxeter arrangement in R^%i\n\n',n);

% Check dimension

  if n < 2
    if verb; fprintf(fout,'\n### data_coxeter: n = %g and must be an integer ≥ 2\n\n',n); end
    V = [];
    return
  end

% Generate

  p = n*(n-1)/2;
  V = zeros(n,p);

  j = 0;

  for i = 1:n-1
    for ii = i+1:n
      j = j+1;
      V(i,j)  =  1;
      V(ii,j) = -1;
    end
  end

return
