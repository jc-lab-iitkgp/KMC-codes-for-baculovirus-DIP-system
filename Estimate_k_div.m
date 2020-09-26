% This code is to estimate the value of the parameter k_div,i.e., 
% Rate of divisions of uninfected cells(CU) and DIP infected cells(CD).
% Without losing generality, in this sub-program we consider only CU.
% Developed by ASHOK DAS
% 09/02/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimated value of parameter k_div = 3.8e-2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

%% Process Parameters
no_cu = 100;  % No. of CU
tot_time = 72;   % Total time in HOURS
time = 0; % Initial time
k_div = 3.8e-2;%        2.31E-2;  % Rate of division of CU and CD;
%% For taking outputs
counter = 1; % counter
delta_t_output = 2; 
t_output = 0; % Initial time
cu_output = no_cu; %Initial CU

%% Process
while time<= tot_time
    event_rate = k_div * no_cu;  % Rate of multiplication of uninfected cells 
    delta_t = 1/ event_rate; % time step
    time = time + delta_t;  % time increament
    no_cu = no_cu +1;  % Cell replication event
    
    % result output
    if time >= (counter* delta_t_output)
        t_output  = [t_output ; time];
        cu_output = [cu_output ; no_cu];
        counter=counter+1;
    end
    % print result
    fprintf('time = %1.4f | no_cu = %1.0f \n', time, no_cu)
end

%%Loading Exp. result
load('CU_density.mat')
figure
plot(t_output, cu_output /cu_output(1),'*-'); % MC
hold on
%plot(time, density,'o-'); % EXP
plot(time(1:4), density(1:4),'o-'); % EXP
legend('MC','Experiment')
title('Cell multiplication: k_{div} = 3.8e-2')
