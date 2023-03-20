function [E_opt] = solver_permuation_ES(N,n,C,C_offset)

E_opt = 99999999;
x_opt = 0;

perms_vec = 1:n;
perms_mat = perms(perms_vec);
for idx_perms = 1:size(perms_mat,1)
    cur_perms_vec = perms_mat(idx_perms,:);
    X = zeros(n,n);
    for idx_row = 1:n
        X(idx_row,cur_perms_vec(idx_row)) = 1;
    end
    
    x = reshape(X,N,1);
    E = x'*C*x + C_offset;
    if(E < E_opt)
        E_opt = E;
        x_opt = x;
    end
    
end

end

