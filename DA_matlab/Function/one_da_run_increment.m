function [Log, x, time] = one_da_run_increment(Problem,Param)
tic;

E_best = 9999999;

% Initialization
x = init_x(Problem,Param);
E_offset = 0;
dE = (2*(1-2*x)).*sum(Problem.Q(:,x==1),2) + diag(Problem.Q);
E_Q_curr = x'*Problem.Q*x;

% Parameter in Tabu
tabu_list = zeros(Problem.N,1);

% Log record
Log.E_Q = [];
Log.flip_idx = [];
Log.E_offset = [];
Log.idx_rand = [];
Log.p_list = [];
Log.step_find_best = [];
Log.E_boost = [];

E_boost = 0;
E_boost_increase_rate = 0.01;

% Start computing
for idx_step = 1:length(Param.temp_sched)
    temp = Param.temp_sched(idx_step);

    p = min(exp(-((dE+E_boost)-E_offset)./temp) ,1);
    accepted = binornd(1,p);
    idx_all = find(accepted);
    if (Param.DA.tabu_tenure)
        idx_tabu = find(tabu_list);
        idx_all = setdiff(idx_all, idx_tabu);
    end
    if( ~isempty(idx_all) )
        idx_rand = randsample(idx_all,1);
        Log.idx_rand(end+1) = idx_rand;
        x(idx_rand) = ~x(idx_rand);
        Log.p_list(end+1) = p(idx_rand);
        E_offset = 0;
        
        E_Q_curr = E_Q_curr + dE(idx_rand);
        if (abs(dE(idx_rand)) < 1e-10 )
            E_boost = E_boost + E_boost_increase_rate;
            E_boost = min(E_boost, Param.DA.E_boost);
        else
            E_boost = 0;
        end
        idx_rand_bar = setdiff(1:Problem.N, idx_rand);
        dE(idx_rand_bar) = dE(idx_rand_bar) - 2*(1-2*x(idx_rand))*(1-2*x(idx_rand_bar)).*Problem.Q(idx_rand,idx_rand_bar)';
        dE(idx_rand) = -dE(idx_rand);
        Log.flip_idx(end+1) = idx_rand;
    else
        E_offset = E_offset + Param.DA.E_offset_increase_rate;
        Log.p_list(end+1) = 0;
        Log.flip_idx(end+1) = -1;
    end

    if (Param.DA.tabu_tenure)
        tabu_list = tabu_list -1; 
        tabu_list(tabu_list<=0) = 0;
        tabu_list(idx_rand) = tabu_list(idx_rand) + Param.DA.tabu_tenure;
    end
    
    if(E_Q_curr + 1e-10 < E_best)
        E_best = E_Q_curr;
        x_best = x;
        Log.step_find_best(end+1) = idx_step;
    end
    
    if(Param.check_Log)
        Log.E_Q(end+1) = E_Q_curr;
        Log.E_offset(end+1) = E_offset;
        Log.E_boost(end+1) = E_boost;
    end
end

x = x_best;

time = toc;
end

