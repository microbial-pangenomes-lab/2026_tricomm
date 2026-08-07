function [R, alpha_C_star, frac] = tricomm_R(P)
% TRICOMM_R  Coexistence ratio R, and the alpha_C that would make R exactly 1.
%
%   [R, alpha_C_star, frac] = tricomm_R(P)
%
% With a = alpha_A/delta, b = alpha_B/delta, c = alpha_C/delta, the no-drug
% model has an interior steady state only when
%
%       a*b*c = a + b + c + 2                                        (*)
%
% R is the ratio of the two sides:
%
%       R = a*b*c / (a + b + c + 2)
%
%   R > 1  -> total biomass grows without bound
%   R < 1  -> the community washes out
%   R = 1  -> a LINE of equilibria (one exactly-zero eigenvalue): the biomass
%             scale is neutral, so the level it settles at is set by the initial
%             condition, not by the parameters.
%
% Note the mu's do not appear. They set the timescale and the transient, but not
% the long-run direction -- see model/PARAMETER_SET_ASSESSMENT.md section 2.
%
% (*) is linear in c, so it can be solved in closed form for alpha_C:
%
%       alpha_C* = delta * (a + b + 2) / (a*b - 1)        (needs a*b > 1)
%
% This is how the wired parameter set was obtained: keep the collaborators'
% fitted mu, alpha_A, alpha_B and delta, and solve for the one remaining alpha.
% alpha_C_star is NaN when a*b <= 1, i.e. no positive alpha_C can balance (*).
%
% frac is the equilibrium biomass composition [LM, PM, PL], from the identity
%   (*)  <=>  1/(1+a) + 1/(1+b) + 1/(1+c) = 1
% whose three terms are exactly the equilibrium fractions. It only sums to 1
% when R == 1; otherwise read it as "the composition the model is heading for".

a = P.alpha_A / P.delta;
b = P.alpha_B / P.delta;
c = P.alpha_C / P.delta;

R = (a * b * c) / (a + b + c + 2);

if a * b > 1
    alpha_C_star = P.delta * (a + b + 2) / (a * b - 1);
else
    alpha_C_star = NaN;
end

frac = [1/(1+a), 1/(1+b), 1/(1+c)];
end
