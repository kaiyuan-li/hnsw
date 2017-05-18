%% initialize
clc; close all;clear all

%% prepare parameters
param.Time = 15e-3;
param.NumOfBeads = 50;
param.Precompression = 0;

% material properties
param.EXX = 195e9*ones(1,param.NumOfBeads);
param.NU = .3*ones(1,param.NumOfBeads);
param.DENS = 7950*ones(1,param.NumOfBeads);
param.G = 9.8;
param.N = 3/2;          % power law coefficient
param.GAMMA = 0;

% deducted parameters
param.DIAMETER = 19e-3*ones(1,param.NumOfBeads);
param.MASS = ones(1,param.NumOfBeads)*4/3*pi.*(param.DIAMETER/2).^3.*param.DENS;
param.K = [inf ...
    1./(3/2/1.414*((1-param.NU(2:param.NumOfBeads).^2)./param.EXX(2:param.NumOfBeads)+(1-param.NU(1:param.NumOfBeads-1).^2)./param.EXX(1:param.NumOfBeads-1))).*...
    sqrt(param.DIAMETER(2:param.NumOfBeads).*param.DIAMETER(1:param.NumOfBeads-1)./(param.DIAMETER(2:param.NumOfBeads)+param.DIAMETER(1:param.NumOfBeads-1)))];
param.DELTA = (param.Precompression./param.K+(1:param.NumOfBeads)*param.MASS(1)*param.G/param.K(2)).^(1/param.N);

param.InitialVelocity = 1;

param.InitialCondition = zeros(1,2*param.NumOfBeads);
param.InitialCondition(2) = param.InitialVelocity;


%% calculate

[t,u] = ode15s(@(t,u)hnswDispFun(t,u,param),[0,param.Time],param.InitialCondition);