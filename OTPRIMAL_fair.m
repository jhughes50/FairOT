function [OutOTDistance,OutPi] = OTPRIMAL_fair(p,q,Theta,omega)
%p = [0.4;0.6];
%q = [0.3;0.2;0.5];
%Theta = [0.1,0.2,0.3;0.2,0.3,0.4];


N = size(p,1);
M = size(q,1);

%Solution
Pi = zeros(N,M);
PiVec = reshape(Pi,[N*M,1]);

lb = zeros(N*M,1);
fun = @(PiVec)Objective(PiVec,Theta,omega,N,M);
func = @(PiVec)NonConstraint(PiVec,p,q);
OutPiVec = fmincon(fun,PiVec,[],[],[],[],lb,[],func);
OutPi = reshape(OutPiVec,[N,M]);
OutOTDistance = OutPiVec'*reshape(Theta,size(OutPiVec)) + omega * sum(log(sum(reshape(OutPiVec,[N,M]),2)+1));

end

function U = Objective(PiVec,Theta,omega,N,M)
U = - PiVec'*reshape(Theta,[size(PiVec,1),1]) - omega * sum(log(sum(reshape(PiVec,[N,M]),2)+1));
end

function [NC,ceq] = NonConstraint(PiVec,p,q)
ceq = [];
NCPi = reshape(PiVec,[size(p,1),size(q,1)]);
NC = [(NCPi*ones(size(q,1),1)-p);(NCPi'*ones(size(p,1),1)-q)];
end