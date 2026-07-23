function ci = nlparci(beta, resid, varargin)
% NLPARCI  Octave shim for MATLAB's nlparci (Statistics Toolbox).
%   CI = nlparci(BETA, RESID, 'Jacobian', J) returns the 95% confidence
%   intervals for the fitted parameters BETA, using the residual vector
%   RESID and the Jacobian J of the model at the solution. CI is
%   [lower, upper], one row per parameter.
%
%   This reproduces MATLAB's linear (Jacobian) approximation:
%       mse    = (r'r) / (n - p)
%       Cov    = mse * inv(J'J)
%       se     = sqrt(diag(Cov))
%       CI     = beta +/- t_{0.975, n-p} * se
%
%   The Student-t quantile is computed from Octave's core betaincinv, so no
%   Statistics package is required (Octave core has no tinv).
%
%   See also: octave_compat/README.md
  J = [];
  for k = 1:2:numel(varargin)
    if strcmpi(varargin{k}, 'jacobian'), J = varargin{k+1}; end
  end
  if isempty(J)
    error('nlparci shim: pass the Jacobian, e.g. nlparci(beta, resid, ''Jacobian'', J).');
  end

  beta  = beta(:);
  resid = resid(:);
  n = numel(resid);
  p = numel(beta);
  dof = max(n - p, 1);

  mse   = (resid' * resid) / dof;
  Sigma = mse * inv(J' * J);
  se    = sqrt(diag(Sigma));
  tcrit = tquant975(dof);
  delta = tcrit .* se;
  ci = [beta - delta, beta + delta];
end

function t = tquant975(v)
% 97.5th percentile of a Student-t distribution with v dof.
% Uses the incomplete-beta inverse: for T ~ t(v) and x = v/(v+t^2),
% the upper-tail prob 1-F(t) = 0.5*I_x(v/2, 1/2). Setting that to 0.025
% gives I_x = 0.05, so x = betaincinv(0.05, v/2, 1/2) and t = sqrt(v(1-x)/x).
  x = betaincinv(0.05, v/2, 0.5);
  t = sqrt(v * (1 - x) / x);
end
