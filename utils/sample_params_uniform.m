function param_sample = sample_params_uniform(param_bounds)
% Generates a single uniformly sampled parameter set.
%
% Arguments:
%   param_bounds - Nx2 matrix of [lower, upper] bounds for each parameter
%
% Returns:
%   param_sample - Nx1 vector of sampled parameter values
% ========================================================================%
lb = param_bounds(:,1);
ub = param_bounds(:,2);

param_sample = lb + (ub - lb) .* rand(numel(lb), 1);
end
