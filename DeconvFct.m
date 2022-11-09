function [gInt,dgInt,TConv,Out,Time] = DeconvFct(f,TExp)
%_______________________________________________________________________________
% [gInt,dgInt,TConv,out,Time] = Deconv(f,Texp)
% Optimization algorithm to perform deconvolution on the experimental data 
% extracted from a TRT to recover a short-term g-function (STgF).The function
% computes the estimated STgF and the estimated convolved
% temperature variations.
%
% Inputs:
%   - f: Incremental temperature function [^oC]
%   - TExp: Experimental temperature variation (T_exp = T_out-T_0) [^oC]
%
% Outputs:
%   - gInt: Interpolated estimated STgF obtained by deconvolution [-]
%   - dgInt: First derivative of the interpolated estimated STgF [-]
%   - TConv: Estimated convolved temperature variation (T_conv=f*ghat) [^oC]
%   - Out: Output variables of the solver [-]
%   - Time: Computing time to converge [s]
%   - global S: Global variable defined to return some parameters for analysis
% 
% Code also validated to deconvolve hydrogeological transfer function under a
% constant pumping test flow rate. In that case, f is the derivative of the
% pumping flow rate and TExp is the experimental drawdown.
%
% Reference: Dion, G., Pasquier, P., & Marcotte, D. (2022). Deconvolution of 
% experimental thermal response test data to recover short-term g-function. 
% Geothermics, 100, 102302. https://doi.org/10.1016/j.geothermics.2021.102302
%
% Author: Gabriel Dion
% Date: 06-2022
% Version: 2.3.0
%_______________________________________________________________________________
%% 0 - Initialization
tic;                                % Start timer
nData = length(TExp);               % Total number of time step
Idall = (1:nData)';                 % Numbered time step signal
DataValidation(f,TExp)              % Inputs validation


%% 1 - Initial guess

% 1.1 - Anchor points specification
Id = unique(round(logspace(log10(1),log10(nData),40)))';
nId = length(Id);

% 1.2 - Initial STgF guess with sequential matrix inversion and fit to an ILS
g0 = ExpIntg0(f,TExp);
g0Id = g0(Id);

% 1.3 - [Temporary] Verification of the initial guess% 
% [TIni] = convFFT(f,g0);
% figure(1001);hold on; plot(TIni); plot(TExp); hold off


%% 2 - Set the solver parameters

% 2.1 - Empty equality constraints
Aeq = [];                           % (empty) Equality constraint
Beq = [];                           % (empty) Equality constraint

% 2.2 - Lower and upper bounds
lb = zeros(nId,1);                  % Lower bound: 0 for positivity constraint
ub = [];

% 2.3 - Linear inequality constraint for positive STgF first derivative
[A,B] = ConstDeriv2(Id);            % C1 and negative second derivative
% [A,B] = ConstDeriv1(nId);         % Positive first derivative constraint
% A = []; B = [];                   % No constraints

% Note: Usually, ConstDeriv2 is best, but if bad results occurs after
% deconvolution, try with ConstDeriv1 and no constraint in last resort.


%% 3 - Weights for the multi-objective function
W.T = 1;                            % Weight for T
W.dg = 10;                          % Weight for dg
W.ddg = 100;                        % Weight for ddg

if nData <= 300
    W.TA = ones(nData,1);
else
    W.TA = [3*ones(300,1);1*ones((nData-300),1)]; % Weight vector for T
end



%% 4 - Set the solver options
opts = optimoptions('fmincon');     % Choose the programming solver
opts.Algorithm = 'interior-point';  % Solving algorithm
opts.MaxFunctionEvaluations = Inf;  % Number of calculation before exiting
opts.MaxIterations = 1000;          % Number of iteration before exiting
opts.OptimalityTolerance = 1e-6;    % Tolerance on the first-order optimality
opts.StepTolerance = 1e-30;         % Termination tolerance
opts.ConstraintTolerance = 0;       % Tolerance on the constraint violation
opts.OutputFcn = {@IterOut};        % Custom function for certain variables
opts.Display = 'final-detailed';    % Display iteration progress symmary
opts.UseParallel = true;            % Use parallel computing
opts.HonorBounds = true;            % Constraints respected at all iteration


%% 5 - Call minimization solver

% Set the equation for the objective function calculation
[gOpt,~,~,Out] = fmincon(@(gIter) ObjFun_ghat(gIter,Id,Idall,f,TExp,W),...
    g0Id,A,B,Aeq,Beq,lb,ub,[],opts);


%% 6 - Interpolation and final convolution

% 5.1 - Interpolation
gInt = pchip(Id,gOpt,Idall);
dgInt = diff([0;gInt]);

% 5.2 - Convolved temperature variations using the estimated STgF (gHat)
TConv = convFFT(f,gInt);

% 5.3 - Output variable
Time = toc;
fprintf("Total optimization time: %2.2f [s]\n",Time)

global S                            % Initiate global variable
S.Id0 = Id;
S.nId0 = nId;
S.g0 = g0;


end


%% Additional functions

function [EGlobal] = ObjFun_ghat(gIter,Id,Idall,f,TExp,W)
%_______________________________________________________________________________
% [EGlobal] = ObjFun_ghat(gIter,Id,Idall,f,TExp,W)
% Computes the multi-objective function used in the deconvolution algorithm.
%
% Inputs:
%   - gIter: Iterated STgF [-]
%   - Id: Position of the nodes n0 to use in interpolation [-]
%   - Idall: Numerized time vector [-]
%   - f: Incremental temperature function [^oC]
%   - TExp: Experimental temperature variation (T_exp=T_out-T_0) [^oC]
%   - W: Structure of weights for the objective function [-]
%
% Output:
%   - EGlobal: Minimized objective function value [-]
%
% Author: Gabriel Dion
% Date: 06-2022
%_______________________________________________________________________________

% 1.1 - Interpolation the STgF to the same length of the time vector
gInt = pchip(Id,gIter,Idall);

% 1.2 - Compute the derivatives of the STgF
dgInt = diff([0;gInt]);             % First derivative
ddgInt = diff([0;gInt],2);          % Second derivative

% 2 - Calculate the convolved temperature variation and partial autocorrelation
TConv = convFFT(f,gInt);            % Convolve temperature

% 3 - Compute each weighted term of the multi-objective function to minimize
E(1) = W.T*rms(W.TA.*(TConv-TExp)); % T
E(2) = W.dg*rms(dgInt);             % dg
E(3) = W.ddg*rms(ddgInt);           % ddg

% 4 - Global value of the multi-objective function to minimize
EGlobal = sum(E);

% 5 - Export the normalized weight of each term of the multi-objective function
global S
S.E = E/EGlobal;
S.EGlobal = EGlobal;
S.TOpt = TConv;
S.ghat = gInt;
end


function [gIni] = ExpIntg0(f,Te)
%_______________________________________________________________________________
% [gIni] = ExpIntg0(f,Te)
% Function that use an exponential integral to generate a first approximation of
% the transfer function used in stationary convolution.
%
% Inputs:
%   - f: Incremental function that describes the test [-] or [degC]
%   - Te: Experimental temperature variation [degC]
%
% Outputs:
%   - gIni: Short-term transfer function initial guess [-]
%
% Author: Gabriel Dion
% Date: 05-2022
%_______________________________________________________________________________

% Vectors of all the same length
if length(f) ~= length(Te)
    error('Input values are not of the same length.')
end

% Initial values to the optimization
Idall = (1:length(Te))';
x0 = [1,1];

% Objective function as the norm l2 of the residual
ObjFun = @(x) rms(convFFT(f,(x(1)*expint(x(2)./Idall)))-Te);

% Performing the minimization of 2 parameters
xOpt = fminsearch(@(x) ObjFun(x),x0);

% Compute the g0 using fitted parameters
gIni = real(xOpt(1)*expint(xOpt(2)./Idall));

% Verification of the presence of complex or not a number in the output
if any(isnan(gIni),'all')
    error('Presence of NaN')
end
end


function [y] = convFFT(f,g)
%_______________________________________________________________________________
% [y]=convFFT(f,g) 
% Performs the non-circular convolution of functions f and g in the spectral 
% domain using the fast Fourier transform and its inverse. The input signals are
% padded with zeros so that the total length of each signal is 2*n-1.
% The convolved signal is then cut to the length of the excitation function f.
%
% Inputs:
%   - f: Excitation function
%   - g: Transfer function of the studied system
%
% Output:
%   - Convolved signal with the same length as the excitation function
%
% Author: Gabriel Dion
% Date: 01-2022
%_______________________________________________________________________________

% Initialization
nData = length(f);
Pad = 2*nData-1;

% Perform non-circular convolution in the spectral domain using padded signals
y = real(ifft(fft(f,Pad).*fft(g,Pad)));

% Cut the signal to the length of the excitation function
y = y(1:nData);

% Verification of the presence of complex or not a number in the output
if any(isnan(y),'all') || any(~isreal(y),'all')
    error('Presence of NaN or complex')
end
end


function [A,B] = ConstDeriv2(Id)
%_______________________________________________________________________________
% [A,B] = NSConstDeriv1(nId)
% Set the linear inequality constraints for the problem Ax <= b used
% in single - deconvolution. The constraints are:
%   1. A positive first derivative constraint (strictly growing function).
%   2. A negative second derivative constraint after an inflexion point.
%
% Input:
%   - Id: Notes position on the transfer function [-]
%
% Outputs:
%   - A: Linear inequality constraint matrix [MxN], where M is the
%       number of inequalities and N is the number of nodes (nId).
%   - B: Vector [Mx1] of constraints to respect.
%
% Author: Gabriel Dion
% Date: 06-2022
%_______________________________________________________________________________

% 0 - Initialization
nId = length(Id);

% 1 - First constraint: positive first derivative

% Vector to set matrix A
e = ones(nId,1);
h = diff([0;Id]);

% Setting the linear inequality equation
A1 = diag(-e)+diag(e(1:(nId-1)),1); % Matrix for first order derivative
A1 = -A1(1:end-1,:)./h(1:end-1);    % Select the [N-1xN] size and correct with h
B1 = zeros(nId-1,1);                % Vector to constrain the first derivative

% 2 - Second constraint: negative second derivative on last section
[~,IdNodes] = min(abs(300-Id));     % Find at around 3h of TRT

A2 = zeros(nId-1-IdNodes,nId);      % Preallocation
for kk = 1:nId-1-IdNodes            % Fill the matrix A2 with second derivative
    c = (Id(IdNodes+kk+1)-Id(IdNodes+kk))/...
        (Id(IdNodes+kk)-Id(IdNodes+kk-1));
    A2(kk,IdNodes+kk-1:IdNodes+kk+1) = [c,-(1+c),1];
end

B2 = zeros(nId-1-IdNodes,1);        % [Mx1] vector with M=N-2

% 3 - Joining the first and second derivative constraint.
A = sparse([A1;A2]);
B = sparse([B1;B2]);

% 4 - Verification
if any(isnan(A),'all') || any(isnan(B),'all')
    error('Presence of NaN')
end
end


function [A,B] = ConstDeriv1(nId)
%_______________________________________________________________________________
% [A,B] = ConstDeriv1(nId)
% Set the linear following inequality constraints for the problem Ax <= b used
% in deconvolution. The constraint is a positive first derivative constraint
% (strictly growing function).
%
% INPUTS:
%   - nId: Total number of points in vector Id [-]
%
% OUPUTS:
%   - A: Linear inequality constraint matrix (MxN), where M is the number of
%       inequalities and N is the number of points (nId).
%   - B: Vector (1xN) of constraints linked to the matrix A.
%
% Author: Gabriel Dion
% Date: 01-2022
%_______________________________________________________________________________

% 1 - Setting the constraint
% Vector to set matrix A
e = ones(nId,1);
% First-order forward finite difference
A = -spdiags([-e,e],0:1,nId-1,nId);
% Zero vector for positive derivative
B = sparse(zeros(nId-1,1));

% 2 - Verification
if any(isnan(A),'all') || any(isnan(B),'all')
    error('Presence of NaN')
end
end


function [stop] = IterOut(~,OptVal,State)
%_______________________________________________________________________________
% Function that save some optimization output variable for each iteration
% in a global variable.
%_______________________________________________________________________________

stop = false;   % Set output so that it never stops the iteration process.

global S

if strcmp(State,'iter')         % Store values at the end of an iteration.
    S.IterTime(OptVal.iteration+1,1) = toc; % Iter time
    S.funcCount(OptVal.iteration+1,1) = OptVal.funccount; % Iter nFunc
    S.feval(OptVal.iteration+1,1) = OptVal.fval; % Iter fval
end
end


function DataValidation(f,TExp)
%_______________________________________________________________________________
% Validate the vectors lengths and their values.
%_______________________________________________________________________________

% Vectors of all the same length
if length(f) ~= length(TExp)
    error('Input values are not of the same length.')
end

% Imaginary value
if any(~isreal(f)) || any(~isreal(TExp))
    error('Input values contain imaginary data point.')
end

% NaN value
if any(isnan(f)) || any(isnan(TExp))
    warning('Input values contain NaN data point.')
    
end

% Infinity value
if any(isinf(f)) || any(isinf(TExp))
    error('Input values contain infiny data point.')
end
end
