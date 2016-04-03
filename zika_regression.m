% MATLAB + JAGS script to perform a hierarchical linear
% regression to estimate initial growth rate (r) for
% Zika virus disease outbreak in Colombia 2015/2016

clear; clc;

% Retrieve observations from CSV file. Skip the first
% row since that row contains column names. After
% retrieval, extract the individual columns as separate
% variables. Note that we transform the case totals
% by taking the natural log.
observations = csvread('dept_totals.csv', 1, 0);
group = observations(:,1);
x     = observations(:,2);
y     = log(observations(:,3));

% Because our data is coded with the initial x-value of 0,
% we have observed intercept values for each group.
% Extact those values so that can be passed to JAGS.
y0 = y(x == 0);
intercept(1:max(group)) = 0;
for i = group(x == 0)
    intercept(i) = y0(find(group(x == 0) == i));
end

% Create the data structure for JAGS.
data = struct( ...
    'group',     group, ...
    'x',         x, ...
    'y',         y, ...
    'intercept', intercept ...
);

% Define initial values for the non-stochastic model
% parameters. Since there are none, this is an empty
% structure.
initial_values(1) = struct;

% Perform the MCMC analysis
[samples, stats, structArray] = matjags( ...
    data, ...
    fullfile(pwd, 'zika_regression.bug'), ...
    initial_values, ...
    'doparallel' ,    0, ...
    'nchains',        1, ...
    'nburnin',        1000, ...
    'nsamples',       10000, ...
    'thin',           1, ...
    'dic',            0, ...
    'monitorparams',  { ...
       'alpha', 'beta', 'beta_mu', 'pvalue' ...
     }, ...
    'savejagsoutput', 1, ...
    'verbosity',      2, ...
    'cleanup',        1, ...
    'rndseed',        0 ...
);

% Use the estimated parameters to calculate predicted values
% and save results in a CSV file.
% for i = 1:length(y)
%     predicted(i) = stats.mean.alpha(group(i)) + ...
%         stats.mean.beta(group(i)) * x(i);
% end
% csvwrite('predicted.csv', [group'; x'; exp(y'); exp(predicted')]');

% Calculate confidence intervals for R0
%
% Heffernan, J. M., Smith, R. J., & Wahl, L. M. (2005).
%   "Perspectives on the basic reproductive ratio." Journal
%   of the Royal Society Interface, 2(4), 281?293.
%   http://doi.org/10.1098/rsif.2005.0042
%
% Zika Virus Disease mean infectious period from Lessler, et
%   al. "Times to Key Events in the Course of Zika Infection
%   and their Implications for Surveillance: A Systematic
%   Review and Pooled Analysis." bioRxiv, 2 Mar 2016.

L = 9.9;

R0 = zeros(0,4);
for i = (group(x == 0))'
    row(1) = i;
    row(2) = 1 + L * stats.ci_low.beta(i);
    row(3) = 1 + L * stats.mean.beta(i);
    row(4) = 1 + L * stats.ci_high.beta(i);
    R0 = [R0; row];
end
% csvwrite('R0.csv', R0);
