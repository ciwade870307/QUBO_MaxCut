function [C,C_offset,P,P_offset,N,E_opt] = load_QAP_file(path,dataname)

disp("Start loading QAP file: "+path+"/qapdata/"+dataname);
data = importdata(path+"/qapdata/"+dataname+".dat");  
sol  = importdata(path+"/qapsoln/"+dataname+".sln");
E_opt = sol(1,2);
data = reshape(data', [],1);
data = data(~isnan(data));
disp("Golden E_opt: "+E_opt)
n = data(1);
N = n*n;    % 2-way 1-hot
F = data(2:N+1); % Flow Matrix
F = reshape(F,n,n);
D = data(N+2:2*N+1); % Flow Matrix
D = reshape(D,n,n);
% Cost
C = kron(F,D);
% Constraint
subP = ones(n,n) - 2*eye(n);
P1 = kron(eye(n),subP); % rowOneHot
P2 = kron(subP,eye(n)); % colOneHot
P = P1+P2;
% Offset
C_offset = 0;
P_offset = 2*n;

end

