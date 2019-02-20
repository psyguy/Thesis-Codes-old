
% Spase Petkoski 03.06.2018, 
% spase.petkoski[at]univ-amu.fr; spase.petkoski[at]gmail.com

% Phase-lags in large scale brain synchronization: Methodological considerations and in-silico analysis

%% wl68N100.mat contains the connectomes (tract weights and lengths) of 100 subjects from HCP; in the article we use subject 50
% It contains variables lnght and wght, for the lengths and weights 
% respectively for the cortical tracts. Each subjects's connectome has 68
% regions. 
%This data has been used for figures 9-15.
%
%% centres.mat contains Cartesian coordinates of the 68 cortical regions from the subject's connectome 
% This data is only used for printing the brain network in Fig. 15.
% 
%% scripts
%% KMcnctm_plv_PlosCB.m and KMfullPlosCB are the scripts used for all the analysis
% The scriptw can be called immediately, and the values for the coupling,
% frequencies, noise, delays distribution, etc. can be changed within it, 
% as described in each script.
% 
%% KMfullPlosCB is used for the analysis shown in figures 1-2 and 5-8
% It computes the time-series of the populations of oscillators (including 
% the case with 2 oscillators) for different frequency and coupling
% distributions. It also calls functions for calculating the phase locking
% values and their statstics.
%
%% KMcnctm_plv_PlosCB is used for the analysis shown in figures 10-15
% It computes the time-series of brain networks consisting of identical 
% phase oscillators connected over the human connectome.
% It also calls functions for calculating the phase locking values and 
% their statstics.
%
%% functions
%
%   plvstatPlosCB --> computes the PLV statistics of the simulated phases
%
%   plvPlosCB --> computes the PLV values
%
%   wtlPlosCB --> computes the weighted lengths for obtaining the mean
%                 inter- and intra-hemispheric lengths
%   
%   KMcnctmHt0PlosCB --> performs the numerical integration for no delays
%
%   KMcnctmHfPlosCB --> performs the numerical integration for cases with 
%                       delays
%
%   alloc_cnctmPlosCB --> allocates the initial values for the simulations
%
%   ind_cnctmPlosCB --> returns indices for the delayed values in the
%                       buffer for each link