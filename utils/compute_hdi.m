function hdi = compute_hdi(data, cred_mass)
% compute_hdi  Compute the Highest Density Interval (HDI) from a sample.
%
%   hdi = compute_hdi(data, cred_mass)
%
%   The HDI is the shortest interval containing the specified probability
%   mass. Unlike equal-tailed intervals, the HDI includes the most probable
%   values and is preferred for skewed or multimodal distributions.
%
%   Inputs:
%     data       - Vector of posterior samples (N x 1)
%     cred_mass  - Credible mass (e.g., 0.89 for 89% HDI)
%
%   Output:
%     hdi        - 1x2 vector [lower, upper] bounds of the HDI

    % Sort samples in ascending order
    sorted_data = sort(data);
    
    % Total number of samples
    n = size(sorted_data, 1);
    
    % Number of samples within the credible interval
    ci_n = floor(cred_mass * n);

    % Width of every possible contiguous interval of length ci_n
    interval_width = sorted_data(ci_n + 1:end) - sorted_data(1:end - ci_n);

    % The HDI is the narrowest such interval
    [~, min_idx] = min(interval_width);

    hdi = [sorted_data(min_idx), sorted_data(min_idx + ci_n)];

end
