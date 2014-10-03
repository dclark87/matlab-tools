%% SIMPLEKALMAN Simple Kalman filter for testing
%
% Author: Daniel Clark, 2013
%
clc;clear all;close all;
%% Load data
% Sample data
T = 5000;
data = .5*randn(1,T)+sin(.004*(1:T));
z = data;                               % init measurement vector
%
zLen = length(z);
%% Init matrices
F = eye(1);                             % state transition matrix
H = eye(1);                             % measurement input matrix
%% Init variables
% Noise paramters:
% PVOLT - initialize your state estimate error
% QVOLT - how much covariance you think the process is varying by (V^2)
% RVOLT - how much covariance in the measurement (sensor noise) (V^2)
PVOLT = 2^2;                            % init P voltage accuracy (covariance) (V^2)
QVOLT = .05^2;                           % covariance of process voltage (V^2)
RVOLT = 1^2;                            % measurement covar error (V^2)
% Covariance matrices
P = diag([PVOLT]);                      % state est. cov. matrix
Q = diag([QVOLT]);                      % process noise cov. matrix
R = diag([RVOLT]);                      % measurement noise cov. matrix
%% Init state and measurement vectors
xhat = zeros(1,zLen);                   % init state vectors
x = zeros(1,zLen);                      % init state vectors
x(1,1) = z(1,1);                        % initial voltage of signal (assume measurement)
% Debugging elements                    % ...
res = zeros(1,zLen);                    % init residuals for plotting
res(:,1) = z(:,1)-H*x(:,1);             % ...
Kres=zeros(1,zLen);                     % Kalman gain * residual
trP = zeros(1,zLen);
trP(1) = P;
%% Iterate through state vector
for k = 2:zLen                          % for each time instant
    % Update state estimate
    xhat(:,k) = F*x(:,k-1);             % update state prediction
    Phat = F*P*F' + Q;                  % estimate P
    % Get measurement z(:,kk)
    K = Phat*H'/(H*Phat*H'+R);          % update Kalman gain
    res(:,k) = z(:,k)-H*xhat(:,k);      % store residuals
    x(:,k) = xhat(:,k)+K*res(:,k);      % state update
    P = (eye(1)-K*H)*Phat;              % update state covariance estimate
    trP(k) = P;
    Kres(:,k) = K*res(:,k);             % store kalman gain'd residuals
    k
end
%% Plots
figure('Name','Plots','Units','Normalized','Position',[0 0 .7 .5]);
subplot(121);plot(z);hold on;plot(x,'r');
title('Measured and filtered signal');legend('Input data (measured)','Filtered data');
subplot(122);plot(trP);title('Trace of estimate covariance (should converge to a number)');
