function [A,a,B,b,x] = data_bdiffmin_dgp_2023

%%
% [A,a,B,b,x] = data_bdiffmin_dgp_2023

  fprintf('Numerical example in section 5.2.9 of the paper');
  fprintf('\n''On the B-differential of the componentwise minimum of two affine vector functions -- The full report''');
  fprintf('\nby J.-P. Dussault, J.Ch. Gilbert and B. Plaqevent-Jourdain [hal-03872711]\n\n');

% Dimension

  n = 3;

% Generate data

  A = eye(n);
  a = zeros(n,1);

  B = [ 2 0 0; 0 2 1; 1 1 2];
  b = zeros(n,1);

  x = zeros(n,1);

return
