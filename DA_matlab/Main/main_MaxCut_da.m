clear; clc; close all
addpath("../Dataset")
addpath("../Function")

figsize = [400,300,700,300];
Param.cvs_en = 0;

%% Load File/  QUBO Formulation/ Penalty weight
instance = 1:1;
dataset = "gset"; % g05/gset/WK2000
path = "../Dataset/MaxCut/"+dataset;

if (Param.cvs_en)
    csvFilename = dataset+instance(1)+"-"+instance(end)+"_da.csv";
    delete(csvFilename)
    csvHeader = ["Instance", "BKS", "Cut_best", "Acc_best", "Cut_avg", "Acc_avg", "Std", "t_avg", "# Opts/ Total runs", "Note"];
    writematrix(csvHeader,csvFilename,'WriteMode','append')
end

for data = instance
dataname = "G"+data+".txt"; % G2.txt
% dataname = "g05_100."+data; % g05_100.0
% dataname = "WK2000_"+data+".rud";
disp("=========================================")
 
[Problem.C, Problem.C_offset, Problem.J, Problem.N, Problem.E_opt] = load_MaxCut_file(path,dataname);
Problem = run_Problem_norm(Problem);

tempEnd_list = [1e-2];
tempStart_list = [1e-1];
E_offset_increase_rate_list = [1];
F_list = [1];
F_decay_list = [1];
E_boost_list = [0];
temp_sched_mode_list = ["exp"];
tabu_tenure_list = [0];

Acc_matrix = zeros(length(tempStart_list),length(tempEnd_list),length(E_offset_increase_rate_list),length(E_boost_list), length(temp_sched_mode_list), length(F_list), length(F_decay_list));

for i = 1:length(tempStart_list)
for j = 1:length(tempEnd_list)
for k = 1:length(E_offset_increase_rate_list)
for l = 1:length(E_boost_list)
for m = 1:length(temp_sched_mode_list)
for n = 1:length(tabu_tenure_list)
for o = 1:length(F_list)
for p = 1:length(F_decay_list)

%% Parameter Setting
% Annealing Parameter
Param.N_run = 1;         % Number of Runs
Param.maxStep = 40000; % 40000; %arg.N^2;    % Number of steps
Param.temp_sched = run_temp_sched(tempStart_list(i),tempEnd_list(j),Param.maxStep,F_list(o),F_decay_list(p),temp_sched_mode_list(m)); % Temperature scheduling
Param.mode.x_init_mode = "rand"; % rand/zeros/ones/spec
Param.mode.x_init = 0;
Param.check_Log = 1;

% DA Parameter
Param.DA.E_offset_increase_rate = E_offset_increase_rate_list(k);
Param.DA.E_boost = E_boost_list(l);    % Used to avoid the ping pong effect when dE = 0
Param.DA.tabu_tenure = tabu_tenure_list(n); % Default: 0

fprintf("-----------------------------------------------\n");
fprintf("tempStart = %f/ tempEnd = %f/ E_offset_increase_rate = %f/ sched = %s/ E_boost = %d/ tabu_tenure = %d/ F = %d/ F_decay = %f\n", ...
        tempStart_list(i), tempEnd_list(j), E_offset_increase_rate_list(k), temp_sched_mode_list(m), E_boost_list(l), tabu_tenure_list(n), F_list(o), F_decay_list(p))

%% Simulation
E_solver_list = [];
time_solver_list = [];
for idx_run = 1:Param.N_run
    [Log, x, time] = one_da_run_increment(Problem,Param);
    % Check result
    Param.E_solver = x'*Problem.C*x;
    E_solver_list(end+1) = Param.E_solver;
    time_solver_list(end+1) = time;
    if(Param.E_solver == Problem.E_opt)
%         fprintf("Find optimal solution!!!\n")
%         break;
    end
end

BKS = abs(Problem.E_opt);
Cut_best = abs(min(E_solver_list));
Acc_best = 100*(min(E_solver_list)/Problem.E_opt);
Cut_avg = abs(mean(E_solver_list));
Acc_avg = 100*mean((E_solver_list)/Problem.E_opt);
Std = abs(std(E_solver_list)/Problem.E_opt*100);
t_avg = mean(time_solver_list);
N_opt = sum(E_solver_list==Problem.E_opt);
Total_run = length(E_solver_list);

fprintf("# of feasible runs: %d\n", Total_run)
fprintf("Avg  E_solver: %f, Sol acc: %f, std: %f\n",Cut_avg, Acc_avg, Std)
fprintf("Best E_solver: %f, Sol acc: %f\n",Cut_best, Acc_best)
fprintf("# of opt sol: %d\n", N_opt)
fprintf("Total time: %f/ Avg time: %f\n", sum(time_solver_list), mean(time_solver_list))

Acc_matrix(i,j,k,l,m,n,o,p) = 100*mean((E_solver_list)/Problem.E_opt);

% Save as csv
if (Param.cvs_en)
% csvHeader = ["Instance", "BKS", "Cut_best", "Acc_best", "Cut_avg", "Acc_avg", "Std", "t_avg", "# Opt run"];
csvrow = [dataname, BKS, Cut_best, Acc_best, Cut_avg, Acc_avg, Std, t_avg, N_opt+"/"+Total_run, ...
                sprintf("tempStart = %f/ tempEnd = %f/ E_offset_increase_rate = %f/ sched = %s/ E_boost = %d/ tabu_tenure = %d/ F = %d/ F_decay = %f", tempStart_list(i), tempEnd_list(j), E_offset_increase_rate_list(k), temp_sched_mode_list(m), E_boost_list(l), tabu_tenure_list(n)), F_list(o), F_decay_list(p)];
writematrix(csvrow,csvFilename,'WriteMode','append')
end

end

end
end
end
end
end
end
end
end

%% Plot Log Result
log_plot_en.temp_sched = 1;
log_plot_en.p_list = 1;
log_plot_en.E_Q = 1;
log_plot_en.flip_idx = 1;
log_plot_en.E_offset = 1;

plot_log_result(Problem, Param, Log, log_plot_en);

% === Energy histogram ===
% histogram(E_solver_list)
% xlabel("Energy")
% ylabel("Frequency")

%% Plot heatmp (tempStart, tempEnd, E_offset_increase_rate_list)
% figure('position',[100,300,1800,400]); hold on; box on; grid on; hold off
% for k = 1:length(E_offset_increase_rate_list)
% subplot(1,length(E_offset_increase_rate_list),k)
% xvalues = string(tempEnd_list);
% yvalues = string(tempStart_list);
% h = heatmap(xvalues,yvalues,Acc_matrix(:,:,k));
% % h = heatmap(xvalues,yvalues,ARPD_matrix, 'ColorLimits',[0.2 2]);
% 
% h.Title = "Acc. @ E_{offset increase rate} = "+E_offset_increase_rate_list(k);
% h.XLabel = 'T_{maxStep-1}';
% h.YLabel = 'T_0';
% end

%% Plot heatmp (tempStart, E_offset_increase_rate_list)
% figure('position',[100,300,1800,400]); hold on; box on; grid on; hold off
% xvalues = string(E_offset_increase_rate_list);
% yvalues = string(tempStart_list);
% h = heatmap(xvalues,yvalues,squeeze(Acc_matrix(:,1,:)));
% % h = heatmap(xvalues,yvalues,ARPD_matrix, 'ColorLimits',[0.2 2]);
% 
% h.XLabel = 'E_{offset_increase_rate_list}';
% h.YLabel = 'T_0';

%% Plot heatmp (F)
% figure('position',[100,300,1000,200]); hold on; box on; grid on; hold off
% xvalues = string(E_offset_increase_rate_list);
% h = heatmap(xvalues,1,squeeze(Acc_matrix(Acc_matrix~=0))');
% h.XLabel = 'E_{offset increase rate}';
