%% Distributed OT, Numerical Experiments
clear
clc
seed = 1;
rng('default')
rng(seed);

MAX_ITER = 250;


% Ex2 N = 5, M = 3
% p3 = rand(10,1);
% p3 = p3./sum(p3);
% q3 = rand(5,1);
% q3 = q3./sum(q3);

N = 5;
M = 2;
p3 = [2;3;4;3;2];
q3 = [4;4];
% Gamma3 = rand(N,M);  % cost parameter on x
% Delta3 = rand(N,M);  % cost parameter on y
Gamma3 = [1,2;3,2;1,4;3,1;2,2];  % cost parameter on x
Delta3 = [2,2;1,2;4,2;2,2;4,4];  % cost parameter on y
Theta3 = Gamma3 + Delta3; % cost parameter on x+y, which is convenient for centralized opt
S3 = ones(N,M);

% weighting consts on fairness
omega1 = 0;  % with fairness
omega3 = 5;  % without fairness


[OTDistance3,Pi3] = OTPRIMALIncomplete_fair(p3,q3,Theta3,S3,omega3); % with fairness
[OTDistance1,Pi1] = OTPRIMALIncomplete_fair(p3,q3,Theta3,S3,omega1); % no fairness
allocation_fair = sum(Pi1,2)
allocation_unfair = sum(Pi3,2)
utility_unfair = OTDistance3
utility_fair = OTDistance1

% iterative algorithm results
[OTDPrimal3,PiPrimal3] = OTRelaxPrimalADMoMPrivateFormal_fair(p3,q3,S3,Gamma3,Delta3,MAX_ITER,omega3); % with fairness
[OTDPrimal1,PiPrimal1] = OTRelaxPrimalADMoMPrivateFormal_fair(p3,q3,S3,Gamma3,Delta3,MAX_ITER,omega1); % no fairness


ReOTDPrimal3 = sqrt(((OTDPrimal3 - OTDistance3*ones(1,MAX_ITER)).^2)); 
ReOTDPrimal1 = sqrt(((OTDPrimal1 - OTDistance1*ones(1,MAX_ITER)).^2)); 

Pi3Re = reshape(Pi3,[1,M*N])';
ResidualPrimal3 = sqrt(sum((PiPrimal3 - Pi3Re*ones(1,MAX_ITER)).^2)); 

Pi1Re = reshape(Pi1,[1,M*N])';
ResidualPrimal1 = sqrt(sum((PiPrimal1 - Pi1Re*ones(1,MAX_ITER)).^2)); 

allocation_alg_final = sum(reshape(PiPrimal3(:,end),[N,M]),2)

%Plots

figure(1)
plot(1:1:MAX_ITER,(OTDistance3*ones(MAX_ITER,1)),'b--',1:1:MAX_ITER,(OTDPrimal3),'b',...
    1:1:MAX_ITER,(OTDistance1*ones(MAX_ITER,1)),'r--',1:1:MAX_ITER,(OTDPrimal1),'r','linewidth',1.5)
hold on
grid on
f1 = legend('Centralized Optimal Solution with Fairness','Distributed Algorithm 1 with Fairness',...
    'Centralized Optimal Solution without Fairness','Distributed Algorithm 1 without Fairness');
f1.FontSize = 12;
xlabel('Iterations')
ylabel('Social Utility under Transport Strategies')
hold off

figure(2)
plot(1:1:MAX_ITER,ResidualPrimal3,'b',1:1:MAX_ITER,ResidualPrimal1,'r-.','linewidth',1.5)
hold on
grid on
f2 = legend('Distributed Algorithm 1 with Fairness','Distributed Algorithm 1 without Fairness');
f2.FontSize = 12;
xlabel('Iterations')
ylabel('\pi Residual')
hold off

figure(3)
plot(1:1:MAX_ITER,ReOTDPrimal3,'b',1:1:MAX_ITER,ReOTDPrimal1,'r-.','linewidth',1.5)
hold on
grid on
f3 = legend('Distributed Algorithm 1 with Fairness','Distributed Algorithm 1 without Fairness');
f3.FontSize = 12;
xlabel('Iterations')
ylabel('Optimal Transport Distance Residual')
hold off

