function [OTDistanceOut,PiOut] = OTRelaxPrimalADMoMPrivateFormal_fair(p,q,S,Gamma,Delta,MAX_ITER,omega)
%p,q are vectors of size $N$, $M$, respectively. 
%p indicates demands, while q indicates supplies.

%S is a matrix of size $N\times M$, which indicates the information of $x$ and $y$. 
%$0$ indicates no infomration between $x$ and $y$, while $1$ indicates information. 

%Gamma and Delta are the utilities of $x$ and $y$, respectively. 

%MAX_ITER indicates the maximum iterations. 

N = size(p,1);
M = size(q,1);

Eta = 1; %We assume eta = 1 for all experiments. 

seed = 1;
rng('default')
rng(seed);

% initialization
Pi1Update = 3*rand(N,M); 
Pi2Update = 3*rand(N,M);
Pi = 3*rand(N,M);
W = rand(N,M);  % updates on alpha in the Alg

for kIterate = 1:MAX_ITER
    %Each x updates
    for nIterate = 1:N
        GammaX = Gamma(nIterate,:);
        pX = p(nIterate);
        SX = S(nIterate,:);
        
        PiX = Pi(nIterate,:);
        WX = W(nIterate,:);
        
        Pi1X = zeros(M,1);   
        fun = @(Pi1X)ObjectivePi1(Pi1X,GammaX,WX,PiX,SX,Eta,omega);
        func = @(Pi1X)NonConstraintPi1(Pi1X,pX,SX);
        OutPi1X = fmincon(fun,Pi1X,[],[],[],[],[],[],func);
        Pi1Update(nIterate,:) = OutPi1X';
    end
    %Each y updates
    for mIterate = 1:M
        DeltaY = Delta(:,mIterate);
        qY = q(mIterate);
        SY = S(:,mIterate);
        
        PiY = Pi(:,mIterate);
        WY = W(:,mIterate);
        
        Pi2Y = zeros(N,1);   
        funm = @(Pi2Y)ObjectivePi2(Pi2Y,DeltaY,WY,PiY,SY,Eta);
        funcm = @(Pi2Y)NonConstraintPi2(Pi2Y,qY,SY);
        OutPi2Y = fmincon(funm,Pi2Y,[],[],[],[],[],[],funcm);
        %VUpdate(mIterate,1) = OutVW(1,1);
        Pi2Update(:,mIterate) = OutPi2Y;
    end
    %Each pair x, y updates
    W = W + (0.5).*Eta.*(Pi1Update-Pi2Update);
    Pi = (0.5).*(Pi1Update+Pi2Update);
    
    PiOut(:,kIterate) = reshape(Pi,[1,N*M]);
    OTDistanceOut(1,kIterate) = reshape((Pi.*S),[1,N*M])*reshape((Gamma+Delta),[1,N*M])'...
        + omega * sum(log(sum(reshape(Pi,[N,M]),2)+1));
end

end


function U = ObjectivePi1(Pi1X,GammaX,WX,PiX,SX,Eta,omega)
U = - (GammaX.*SX)*Pi1X + (WX.*SX)*Pi1X + (Eta./2).*(Pi1X'.*SX - PiX.*SX)*(Pi1X'.*SX - PiX.*SX)'...
    - omega*log(SX*Pi1X+1);  % with fairness here
end

function U = ObjectivePi2(Pi2Y,DeltaY,WY,PiY,SY,Eta)
U = - (DeltaY.*SY)'*Pi2Y - (WY.*SY)'*Pi2Y + (Eta./2).*(Pi2Y.*SY - PiY.*SY)'*(Pi2Y.*SY - PiY.*SY);
end

function [NC,ceq] = NonConstraintPi1(Pi1X,pX,SX)
NC1 = ones(size(Pi1X))'*(Pi1X.*SX')  - pX ;
NC = [NC1;(-Pi1X.*SX')];
ceq = [];
end

function [NC,ceq] = NonConstraintPi2(Pi2Y,qY,SY)
NC2 =  ones(size(Pi2Y))'*(Pi2Y.*SY)  - qY ;
NC = [NC2;(-Pi2Y.*SY)];
ceq = [];
end