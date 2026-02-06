function [y, x] = one_state_simulation_with_noise(parameters, r, EC, EC_value, SB)
% Adapted from LMSE v1.2 by Scott Albert
% Link: https://storage.googleapis.com/wzukusers/user-31382847/documents/5a60fe4a59003gmGKILL/LMSE.v1.2.zip
%
% Author: Josh Liddy
% Email: jliddy@umass.edu
% Institution: UMass Amherst
% Date: February 18, 2025
% Version 1.0
%
%
% Summary:
%    This function simulates a single-rate state space model.
%
% Arguments:
%    parameters: the two-state model parameters
%    r: the perturbation on each trial
%    EC: an array that indicates if a trial is an error-clamp trial
%        If the n-th entry is non-zero, this indicates that trial n is an
%        error-clamp trial
%        If the n-th entry is zero, this indicates that trial n is not an
%        error-clamp trial
%    EC_value: an array that indicates the value of the clamped error on
%        each error-clamp trial. 
%    SB: an array that indiciates if a trial is followed by a set break
%        If the n-th entry is non-zero, this indicates that trial n is
%        followed by a set break
%        If the n-th entry is zero, this indicates that trial n is not
%        followed by a set break
%
% Returns:
%    y: motor output
%    x: state estimate
% ======================================================================= %

% Parameters
a = parameters(1);     % Retention
b = parameters(2);     % Error sensitivity
x1 = parameters(3);    % Initial state
d = parameters(4);     % Decay exponent
sigma = parameters(5); % Measurement noise SD

% Number of trials
N = length(r);

% Preallocate state and output arrays and initialize.
x = zeros(N,1);
y = zeros(N,1);
x(1) = x1;
y(1) = x(1);

% Simulate trials
for n = 2:N

    % Determie error based on error clamp
    if EC(n-1) == 0
        e = r(n-1) - y(n-1);
    else
        e = EC_value(n-1);
    end

    % Determine retention based on set break
    if SB(n-1) == 0
        A = a;
    else
        A = a^d;
    end

    % Update state
    x(n) = A * x(n-1) + b * e;

    % Generate output
    y(n) = x(n);
end

% Add measurement noise to y
y = y + normrnd(0,sigma,size(y));

end