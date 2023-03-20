function [w] = run_penalty_weight(C,G,N)
%% Refer to "Penalty Weights in QUBO formulations: Permutation Problems"

w = {};
C_triu = 2*triu(C,1) + diag(diag(C)); % Symmetric matrix to upper triangle matirx 
G_triu = 2*triu(G,1) + diag(diag(G)); % Symmetric matrix to upper triangle matirx 
%% UB
w.ub = sum(abs(C_triu),'all');  
%% MQC
% w.mqc = max(C_triu, [], 'all');  
%% VLM
Wc = [];
for i = 1:N
    % min
    idx_min = C_triu(i,:)<0;
    idx_min(i) = 0;
    Wc(end+1) = -C(i,i)-sum(C_triu(i,idx_min));
    % max
    idx_max = C_triu(i,:)>0;
    idx_max(i) = 0;
    Wc(end+1) = C(i,i)+sum(C_triu(i,idx_max));
end
w.vlm = max(Wc);
%% MOMC
Wg = [];
for i = 1:N
    % min
    idx_min = G_triu(i,:)<0;
    idx_min(i) = 0;
    Wg(end+1) = -G(i,i)-sum(G_triu(i,idx_min));
    % max
    idx_max = G_triu(i,:)>0;
    idx_max(i) = 0;
    Wg(end+1) = G(i,i)+sum(G_triu(i,idx_max));
end
gamma = min(Wg(Wg>0));
w.momc = max(1, w.vlm/gamma);
%% MOC
W_abs_ratio = abs(Wc./Wg);
w.moc = max(1,max(W_abs_ratio((Wg>0))));

end

