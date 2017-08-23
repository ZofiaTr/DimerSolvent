%% Compute variance from ACF


clear all;
clc;
clear all hidden;
close all;


%fileBase='C:\Users\trstanova\Documents\gForge\Simulations\RESULTS\Variance Dimer Solvent Second Order Kmin 0-2.7 Kmax 3';
%fileData='BigEps';

%cd(fileBase);



%% load data
AUTOCORRELATION_FUNCTION_1;  % A1
VARIANCE_IN_NUMBER_OF_REPLICAS_1; % Variance1


AUTOCORRELATION_FUNCTION_2;
VARIANCE_IN_NUMBER_OF_REPLICAS_2;


AUTOCORRELATION_FUNCTION_3;
VARIANCE_IN_NUMBER_OF_REPLICAS_3;


%VARIANCE_IR; %only a number

%%


%
%

ACF1=A1-AAverEnd1^2;
ACF2=A2-AAverEnd2^2;
ACF3=A3-AAverEnd3^2;

Variance1B =  2*dt*trapz(ACF1);
Variance2B =  2*dt*trapz(ACF2);
Variance3B =  2*dt*trapz(ACF3);

%         Variance1B =  2*dt*trapz(A1-AAver1.*AAver1(1));
%         Variance2B =  2*dt*trapz(A2-AAver2.*AAver2(1));
%         Variance3B =  2*dt*trapz(A3-AAver3.*AAver3(1));
%
%

fprintf('Variance 1 = %f \n', Variance1B);
fprintf('Variance 2 = %f \n', Variance2B);
fprintf('Variance 3 = %f \n', Variance3B);


fprintf('Variance 1 from C++ = %f \n', Variance1);
fprintf('Variance 2 from C++ = %f \n', Variance2);
fprintf('Variance 3 from C++ = %f \n', Variance3);


figure(1)
plot(ACF1, '-r')
figure(2)
plot(ACF2, '-b')
figure(3)
plot(ACF3, '-g')
