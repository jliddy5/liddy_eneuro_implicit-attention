function [parameters, mse] = LMSE(y, r, EC, EC_value, SB, param_init, param_bound)
% Adapted from LMSE v1.2 by Scott Albert
% Link: https://storage.googleapis.com/wzukusers/user-31382847/documents/5a60fe4a59003gmGKILL/LMSE.v1.2.zip
%
% Author: Josh Liddy
% Email: jliddy@umass.edu
% Institution: UMass Amherst
% Date: February 18, 2025
% Version 1.0
%
% Summary:
%    This function minimizes the squared error between observed data and
%    predictions from a single-rate state-space model. Equivalent to
%    maximum likelihood estimation under assumptions of noise in movement
%    production but not planning.
%
% Arguments:
%    y: Observed motor output
%    r: Rotation sequence
%    EC: Binary array indicating error-clamp trials
%    EC_value: Error value during error-clamp trials
%    SB: Binary array indicating set breaks
%    param_init: Initial guess for model parameters [a, b, x1, d]
%    param_bound: Lower and upper bounds for parameters
%
% Returns:
%    parameters: Estimated parameters [a, b, x1, d]
%    mse: Mean squared error of the fitted model
% ======================================================================= %

    function MSE = cost_fun(params)

        % Check if parameters exceed bounds and penalize violations.
        if any(params < lb) || any(params > ub)
            MSE = Inf;
            return;
        end

        % Fixed initial state x(1) = 0
        all_params = [params(1:2); 0; params(3)];

        % Simulate model
        yPred = one_state_simulation_without_noise(all_params, r, EC, EC_value, SB);

        % Compute the mean squared error
        MSE = mean((y-yPred).^2,"omitmissing");
    end

% Parameter bounds
lb = param_bound(:,1);
ub = param_bound(:,2);

% Settings for optimization
options = optimset('Display', 'off');

% Perform nonlinear optimization with fmincon. Add in fixed x(1) = 0.
parameters = fmincon(@cost_fun, param_init, [], [], [], [], lb, ub, [], options);
parameters = [parameters(1:2); 0; parameters(3)];

% Compute mean-squared error
yPred = one_state_simulation_without_noise(parameters, r, EC, EC_value, SB);
mse = mean((y - yPred).^2,"omitmissing");

end