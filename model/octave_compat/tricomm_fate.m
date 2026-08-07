function [s, code] = tricomm_fate(R, tol)
% TRICOMM_FATE  One-word long-run verdict from the coexistence ratio R.
%
%   [s, code] = tricomm_fate(R)          tol defaults to 1e-4
%
% code is the short form used in sweep tables: 'grow', 'washout', 'neutral'.
% s spells it out. R = 1 is neutral rather than stable on purpose: the model has
% a LINE of equilibria there, so the level reached depends on the initial
% condition (see tricomm_R and PARAMETER_SET_ASSESSMENT.md section 7).
%
% WHY tol = 1e-4. Near R = 1 the total biomass drifts at very close to
%   8 * (R - 1)  decades per 1000 h
% (measured over alpha_C = 0.055 .. 0.071 on calc_par; the coefficient is
% constant to ~7% across that range). So |R - 1| = 1e-4 moves the community by
% 0.2% over an 1100 h run and ~10% over 50 000 h -- far below anything an
% experiment could resolve, and worth calling neutral. Pass a smaller tol if you
% specifically want to see how exactly a set sits on the surface; the wired
% solved.par has R - 1 = 5e-6, i.e. 0.4% over 50 000 h.

if nargin < 2, tol = 1e-4; end

if ~isfinite(R)
    code = 'invalid';
    s    = 'not computable';
elseif R > 1 + tol
    code = 'grow';
    s    = 'unbounded growth';
elseif R < 1 - tol
    code = 'washout';
    s    = 'washout';
else
    code = 'neutral';
    s    = 'neutral (line of equilibria; level set by the initial condition)';
end
end
