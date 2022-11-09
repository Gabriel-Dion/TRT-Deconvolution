% DECONVOLUTION OF EXPERIMENTAL THERMAL RESPONSE TEST DATA TO RECOVER 
% SHORT-TERM g-FUNCTION - Algorithm
%________________________________________________________
% DESCRIPTION OF THE SCRIPT
%
% This script, in conjunction with the associated functions and simulation data,
% allows to perform a deconvolution algorithm to the entire data set of a TRT 
% in order to extract a short-term g-function (STgF), its first derivative and 
% the convolved temperature reconstruction.
%________________________________________________________
% REFERENCE
%
% Refere to this work as:
%
% Dion, G., Pasquier, P., & Marcotte, D. (2021). DECONVOLUTION OF EXPERIMENTAL
% THERMAL RESPONSE TEST DATA TO RECOVER SHORT-TERM g-FUNCTION. Geothermics, 99.
%
%________________________________________________________
% Author: Gabriel Dion
%________________________________________________________
%% 1 - Import synthetic Data
clearvars; clc;

% 1.1 - Add paths containing and datasets
addpath(strcat(pwd,'/Datasets'));

% 1.2 - Import the TRT data set variables [either 1 or 2]
nTest = 2;
[tData,nData,gRef,V,Q,f,Tin,Tout,T] = ImportNum(nTest);


%% 2 - Deconvolution algorithm

% 2.1 - Deconvolution
[g,dg,TConv,Out,Time] = DeconvFct(f,T);

% 2.1 - Printing some outputs
rmseT = rms(TConv-T);
rmseg = rms(g-gRef);

fprintf('Temperature RMSE: %2.2f [degC]\n',rmseT)
fprintf('STgF RMSE: %2.2f [degC]\n',rmseg)


%% 3 - Figure

figure(1)

subplot(1,2,1)

yyaxis left;
hold on
plot(tData,g,'LineWidth',1.2,'Color','b')
plot(tData,gRef,'LineWidth',1,'Color','k','LineStyle','--')
hold off

ax = gca; ax.FontName = 'Times'; ax.FontSize = 9;
ax.YMinorTick = 'on'; ax.XMinorTick = 'on';
ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k'; ax.XScale = 'log';

yyaxis right;
hold on
plot(tData,dg,'LineWidth',1.2,'Color','b')
plot(tData,diff([0;gRef]),'LineWidth',1,'Color','k','LineStyle','--')
hold off

xlabel('Time (d)','interpreter','latex');
yyaxis left; ylabel('$\hat{g}$ [-]','interpreter','latex')
yyaxis right; ylabel('$\hat{g}''$ [-]','interpreter','latex')
legend('$\hat{g}$','$g$','$\hat{g}$''','$g$''',...
    'location','best','fontsize',9,'fontname','Times','interpreter','latex')
box on; ax = gca; ax.FontName = 'Times'; ax.FontSize = 9;
ax.YMinorTick = 'on'; ax.XMinorTick = 'on';
ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k'; ax.XScale = 'log';

subplot(1,2,2)

hold on
plot(tData,TConv,'LineWidth',1.2,'Color','b')
plot(tData,T,'LineWidth',1,'Color','k','LineStyle','--')
hold off

xlabel('Time (d)','interpreter','latex');
ylabel('$T_{out}-T_0$ [$^oC$]','interpreter','latex')
legend('Convolved','Experimental',...
    'location','best','fontsize',9,'fontname','Times','interpreter','latex')
box on; ax = gca; ax.FontName = 'Times'; ax.FontSize = 9;
ax.YMinorTick = 'on'; ax.XMinorTick = 'on';
