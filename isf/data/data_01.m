function [V] = signed_feasibility_data01(n,m)

%%
% [V] = signed_feasibility_data01(n,m)
%
% Provides random vectors for signed_feasibility_main. The ouput matrix V is mxn.

  rng(0);		% set the seed of the random generator
  V = rand(n,m)-0.5;	% generate the vectors into the columns of the matrix V

return
