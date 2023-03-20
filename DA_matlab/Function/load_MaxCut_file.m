function [C,C_offset,J,N,E_opt] = load_MaxCut_file(path,dataname)

disp("Start loading MaxCut file: "+path+"/"+dataname);
sol = importdata(path+"/solution.txt");
E_opt = -sol.data(sol.textdata == dataname);
disp("Golden E_opt: "+(E_opt))
data = importdata(path+"/"+dataname);  
data = reshape(data', [],1);
data = data(~isnan(data));
N = data(1);
N_edge = data(2);
data = reshape(data(3:end),3,N_edge)';
J = zeros(N,N);
C = zeros(N,N);
C_offset = 0;
for i = 1:N_edge
    J( data(i,1), data(i,2) ) = data(i,3);
%     C( data(i,1), data(i,2) ) = C( data(i,1), data(i,2) ) + data(i,3);
%     C( data(i,2), data(i,1) ) = C( data(i,2), data(i,1) ) + data(i,3);
%     C( data(i,1), data(i,1) ) = C( data(i,1), data(i,1) ) - data(i,3);
%     C( data(i,2), data(i,2) ) = C( data(i,2), data(i,2) ) - data(i,3);
end
C = J + J';
C = C - diag(sum(J,1)'+sum(J,2));

% [C, C_offset] = ising2qubo(J,zeros(N,1));

end

