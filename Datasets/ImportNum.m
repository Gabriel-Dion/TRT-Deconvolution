function [tData,nData,gRef,V,Q,f,Tin,Tout,T] = ImportNum(nTest)
%_______________________________________________________________________________
% [tData,nData,g,V,Q,f,T0,EWT,LWT,dT,Stamp] = ImportNum(nTest)
% Function that import TRT simulation test cases based on the numerical model 
% of Beaudry et al. 2019.
%
% Input:
%   - nTest[1x1]: Number of the test to import [-] 
%
% Outputs:
%   - tData: Time step vector [d]
%   - nData: Number of time step in time array [-]
%   - g: Refence transfer function obtain with numerical model [-]
%   - V: Circulation flow rate [m^3/s]
%   - Q: Heating power [W]
%   - f: Incremental temperature function [^oC]
%   - Tin: Inlet water temperature [^oC]
%   - Tout: Outlet water temperature [^oC]
%   - T: Temperature variation Tout - T0 [^oC]
%
% REFERENCE:
%
% Beaudry, Gabrielle, Philippe Pasquier, and Denis Marcotte. ‘The Impact of 
% Rock Fracturing and Pump Intake Location on the Thermal Recovery of a Standing
% Column Well: Model Development, Experimental Validation, and Numerical 
% Analysis’. Science and Technology for the Built Environment 25, no. 8 
% (14 September 2019): 1052–68. https://doi.org/10.1080/23744731.2019.1648133.
%
% Author: Gabriel Dion
% Date: 04-2021
%_______________________________________________________________________________

% Initial parameters (numerical model's input)
tData = linspace(0,7,7*24*60)';             % Minutes over 7 days [d]
nData = length(tData);                      % Number of data [-]

T0 = 11;                                    % Undisturbed ground temperature [^oC]
V = (100/60000);                            % 100 L/min Flow rate [m^3/s]

% Reference STgF
GFUNC = csvread('gRef.csv',5,0);            % Reference STgF
gRef = pchip(GFUNC(:,1),GFUNC(:,2),tData);  % Interpolated STgF [-]

% Heating power
HEATFLUX = csvread(strcat('Q',num2str(nTest),'.csv'));
ERR.HEAT = csvread('NoiseHeat.csv');
Q = pchip(HEATFLUX(:,1),HEATFLUX(:,2),tData)+...
    pchip(ERR.HEAT(:,1),ERR.HEAT(:,3),tData)/2;

% Temperature
TEMP = csvread(strcat('T',num2str(nTest),'Simul.csv'),5,0);
ERR.TEMP = csvread('NoiseTemp.csv');

Tin = pchip(TEMP(:,1),TEMP(:,3),tData)+...
    pchip(ERR.TEMP(:,1),ERR.TEMP(:,3),tData)/2;
Tout = pchip(TEMP(:,1),TEMP(:,2),tData)+...
    pchip(ERR.TEMP(:,1),ERR.TEMP(:,3),tData)/2;

% Temperature variation
T = Tout-T0;

% Incremental temperature function
Cp_w = 4183;                                % Specific water capacity [J/kgK]
rho_w = 999;                                % Water density [kg/m^3]
f = diff([0;(Q)/(V*Cp_w*rho_w)]);

end

