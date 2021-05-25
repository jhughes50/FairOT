%% Distributed OT, Numerical Experiments
clear
clc
seed = 1;
rng('default')
rng(seed);

Pi_initial_primal = ones(3,2);
W_initial_primal = zeros(3,2);


MAX_ITER = 250;

omega = 4; % fairness coefficient
%Ex1
p1 = [5;10;2];
q1 = [5;10];
Gamma1 = [2,0;
    3,4;
    0,4];
Delta1 = [5,0;
    4,4;
    0,3];
Theta1 = Gamma1 + Delta1;
S1 = [1,0;
    1,1
    0,1];
[OTDistance1,Pi1] = OTPRIMALIncomplete_fair(p1,q1,Theta1,S1,omega);
[OTDPrimal1,PiPrimal1,WPrimal1] = OTRelaxPrimalADMoMPrivateFormalOnline_fair(p1,q1,S1,Gamma1,Delta1,MAX_ITER,Pi_initial_primal,W_initial_primal,omega); % online
[OTDPrimal1_off,PiPrimal1_off] = OTRelaxPrimalADMoMPrivateFormal_fair(p1,q1,S1,Gamma1,Delta1,MAX_ITER,omega); % not online; random initialization in the function

Pi1Re = reshape(Pi1,[1,6])';

ReOTDPrimal1 = sqrt(((OTDPrimal1 - OTDistance1*ones(1,MAX_ITER)).^2)); 
ResidualPrimal1 = sqrt(sum((PiPrimal1 - Pi1Re*ones(1,MAX_ITER)).^2));

ReOTDPrimal1_off = sqrt(((OTDPrimal1_off - OTDistance1*ones(1,MAX_ITER)).^2)); 
ResidualPrimal1_off = sqrt(sum((PiPrimal1_off - Pi1Re*ones(1,MAX_ITER)).^2));


% Pi_initial_primal = reshape(PiPrimal1(:,MAX_ITER),[3,2]);
% W_initial_primal = reshape(WPrimal1(:,MAX_ITER),[3,2]);

Pi_initial_primal = [reshape(PiPrimal1(:,MAX_ITER),[3,2]),ones(3,1)];
W_initial_primal = [reshape(WPrimal1(:,MAX_ITER),[3,2]),zeros(3,1)];

%Ex2
% p2 = [8;6;4];
% q2 = [5;8];
% Gamma2 = [4,2;
%     3,2;
%     0,3];
% Delta2 = [6,4;
%     3,5;
%     0,6];
% Theta2 = Gamma2 + Delta2;
% S2 = [1,1;
%     1,1
%     0,1];
% [OTDistance2,Pi2] = OTPRIMALIncomplete_fair(p2,q2,Theta2,S2,omega);
% [OTDPrimal2,PiPrimal2,WPrimal2] = OTRelaxPrimalADMoMPrivateFormalOnline_fair(p2,q2,S2,Gamma2,Delta2,MAX_ITER,Pi_initial_primal,W_initial_primal,omega);
% 
% 
% Pi2Re = reshape(Pi2,[1,6])';
% 
% ReOTDPrimal2 = sqrt(((OTDPrimal2 - OTDistance2*ones(1,MAX_ITER)).^2)); 
% 
% ResidualPrimal2 = sqrt(sum((PiPrimal2 - Pi2Re*ones(1,MAX_ITER)).^2)); 
% 
% 
% Pi_initial_primal = [reshape(PiPrimal2(:,MAX_ITER),[3,2]),ones(3,1)];
% W_initial_primal = [reshape(WPrimal2(:,MAX_ITER),[3,2]),zeros(3,1)];


%Ex3
p3 = [8;6;14];
q3 = [7;8;5];
Gamma3 = [4,2,0;
    3,2,5;
    0,3,3];
Delta3 = [6,4,0;
    3,5,4;
    0,6,2];
Theta3 = Gamma3 + Delta3;
S3 = [1,1,0;
    1,1,1
    0,1,1];

[OTDistance3,Pi3] = OTPRIMALIncomplete_fair(p3,q3,Theta3,S3,omega);
[OTDPrimal3,PiPrimal3,WPrimal3] = OTRelaxPrimalADMoMPrivateFormalOnline_fair(p3,q3,S3,Gamma3,Delta3,MAX_ITER,Pi_initial_primal,W_initial_primal,omega);
[OTDPrimal3_off,PiPrimal3_off] = OTRelaxPrimalADMoMPrivateFormal_fair(p3,q3,S3,Gamma3,Delta3,MAX_ITER,omega); % not online; random initialization in the function

Pi3Re = reshape(Pi3,[1,9])';

ReOTDPrimal3 = sqrt(((OTDPrimal3 - OTDistance3*ones(1,MAX_ITER)).^2)); 
ResidualPrimal3 = sqrt(sum((PiPrimal3 - Pi3Re*ones(1,MAX_ITER)).^2)); 

ReOTDPrimal3_off = sqrt(((OTDPrimal3_off - OTDistance3*ones(1,MAX_ITER)).^2)); 
ResidualPrimal3_off = sqrt(sum((PiPrimal3_off - Pi3Re*ones(1,MAX_ITER)).^2)); 

Pi_initial_primal = reshape(PiPrimal3(:,MAX_ITER),[3,3]);
W_initial_primal = reshape(WPrimal3(:,MAX_ITER),[3,3]);



%Ex4
p4 = [8;10];
q4 = [8;5;7];
Gamma4 = [3,2,5;
    0,3,3];
Delta4 = [3,5,4;
    0,6,2];
Theta4 = Gamma4 + Delta4;
S4 = [1,1,1
    0,1,1];

[OTDistance4,Pi4] = OTPRIMALIncomplete_fair(p4,q4,Theta4,S4,omega);
[OTDPrimal4,PiPrimal4,WPrimal4] = OTRelaxPrimalADMoMPrivateFormalOnline_fair(p4,q4,S4,Gamma4,Delta4,MAX_ITER,Pi_initial_primal(2:3,:),W_initial_primal(2:3,:),omega);
[OTDPrimal4_off,PiPrimal4_off] = OTRelaxPrimalADMoMPrivateFormal_fair(p4,q4,S4,Gamma4,Delta4,MAX_ITER,omega); % not online; random initialization in the function

Pi4Re = reshape(Pi4,[1,6])';

ReOTDPrimal4 = sqrt(((OTDPrimal4 - OTDistance4*ones(1,MAX_ITER)).^2)); 
ResidualPrimal4 = sqrt(sum((PiPrimal4 - Pi4Re*ones(1,MAX_ITER)).^2));

ReOTDPrimal4_off = sqrt(((OTDPrimal4_off - OTDistance4*ones(1,MAX_ITER)).^2)); 
ResidualPrimal4_off = sqrt(sum((PiPrimal4_off - Pi4Re*ones(1,MAX_ITER)).^2));

%%

%Plot
figure(1)
plot(1:1:(3*MAX_ITER),[OTDistance1*ones(MAX_ITER,1);OTDistance3*ones(MAX_ITER,1);OTDistance4*ones(MAX_ITER,1)],'k',...
    1:1:(3*MAX_ITER),[OTDPrimal1,OTDPrimal3,OTDPrimal4],'b--',1:1:(3*MAX_ITER),[OTDPrimal1_off,OTDPrimal3_off,OTDPrimal4_off],'r-.','linewidth',1.5)
hold on
grid on
legend('Centralized Optimal Solution with Fairness','Online Distributed Algorithm with Fairness','Distributed Algorithm 1 with Fairness')
xlabel('Iterations')
ylabel('Social Utility under Transport Strategies')
hold off

figure(2)
plot(1:1:(3*MAX_ITER),[ResidualPrimal1,ResidualPrimal3,ResidualPrimal4],'b',...
    1:1:(3*MAX_ITER),[ResidualPrimal1_off,ResidualPrimal3_off,ResidualPrimal4_off],'r-.','linewidth',1.5)
hold on
grid on
legend('Online Distributed Algorithm with Fairness','Distributed Algorithm 1 with Fairness')
xlabel('Iterations')
ylabel('\pi Residual')
hold off

figure(3)
plot(1:1:(3*MAX_ITER),[ReOTDPrimal1,ReOTDPrimal3,ReOTDPrimal4],'b',...
    1:1:(3*MAX_ITER),[ReOTDPrimal1_off,ReOTDPrimal3_off,ReOTDPrimal4_off],'r-.','linewidth',1.5)
hold on
grid on
legend('Online Distributed Algorithm with Fairness','Distributed Algorithm 1 with Fairness')
xlabel('Iterations')
ylabel('Optimal Transport Distance Residual')
hold off


