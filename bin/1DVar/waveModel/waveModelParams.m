function [g,alpha,beta,nu,gammaType]=waveModelParams();
%
% common params for NL-TL-AD model
%

g=9.8;
alpha=1;
beta=0.1;  % roller parameter
nu=.5;

% switch for breaker model.  Use '2001' Battjes and Stive (1985) (as used by
% Ruessink et al. 2001), or '2003' for Ruessink et al. (2003).
gammaType=2003;
