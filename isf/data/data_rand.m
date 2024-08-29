function [V] = signed_feasibility_data_rand(n,m,r,fout,verb)

%%
% [V] = signed_feasibility_data_rand(n,m,r,fout,verb)
%
% Provides random vectors for signed_feasibility_main. The ouput matrix
% V is mxn of rank r.

  if verb; fprintf(fout,'\nRandom problem\n\n'); end

% Check dimension

  if (r < 1) | ( r > min(n,m) ) | (fix(r)-r)
    if verb; fprintf(fout,'\n###Random problem: r=%i must be an integer in [1:min(%i,%i)]\n\n',p,n,m); end
    V = [];
    return
  end

% Generate

  rng(0);				% set the seed of the random generator

  if r < n
    V = data_01(n,r);		% r random vectors used to determine a range space of dimension r
    B = orth(V);		% orthogonal basis of R(V)
    V = data_01(n,m);		% random vectors (args (n,m) => V is nxm)
    V = B*((B'*B)\(B'*V));	% project V on R(B)
  else
    V = data_01(n,m);		% random vectors (args (n,m) => V is nxm)
  end

return
