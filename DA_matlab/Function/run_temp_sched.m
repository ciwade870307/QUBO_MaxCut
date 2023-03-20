function [temp_sched] = run_temp_sched(tempStart,tempEnd,maxStep,F,F_decay,mode)

% ========================================================================
if mode == "exp"
    tempDecay = (tempEnd/tempStart)^(1/(maxStep-1));  % Temperature Decay
    temp_sched = tempStart*(tempDecay.^(0:maxStep-1));
    tmp = temp_sched(1:F:end); % Field-Cycling, Default: 1 e.g, no Field-Cycling
    % arg.temp_sched = repmat(tmp, 1, arg.F);
%     F_decay = 0.7;
    for f = 1:F
        temp_sched(1+(f-1)*maxStep/F:(f)*maxStep/F) = (tmp-tempEnd)*(F_decay)^(f-1)+tempEnd;
    end
    
    % N_loop = 1;
    % tmp = arg.temp_sched(1:N_loop:end);
    % for idx_loop = 1:N_loop
    %     arg.temp_sched(idx_loop:N_loop:end) = tmp;
    % end

% ========================================================================
elseif mode == "constExp"   
    R_constExp = 0.2;
    conStep = R_constExp*maxStep;       % # Constant steps
    linaerStep = (1-R_constExp)*maxStep;   % # Exponential steps
    temp_sched(1:conStep) = tempStart;
    tempDecay = (tempEnd/tempStart)^(1/(linaerStep-1));  % Temperature Decay
    temp_sched(conStep+1:maxStep) = tempStart*(tempDecay.^(0:linaerStep-1));

% ========================================================================
elseif mode == "constLinear"   
    R_constExp = 0.2;
    conStep = R_constExp*maxStep;       % # Constant steps
    linaerStep = (1-R_constExp)*maxStep;   % # Linear steps
    temp_sched(1:conStep) = tempStart;
    temp_sched(conStep+1:maxStep) = linspace(tempStart,tempEnd,linaerStep);
    
% ========================================================================
elseif mode == "log"    
    denom = (1:maxStep)+1000;
    temp_sched = 1./log(denom);
    
    
% ========================================================================
elseif mode == "linear"
    temp_sched = linspace(tempStart,tempEnd,maxStep);
    tmp = temp_sched(1:F:end); % Field-Cycling, Default: 1 e.g, no Field-Cycling
%     F_decay = 1;
    for f = 1:F
        temp_sched(1+(f-1)*maxStep/F:(f)*maxStep/F) = (tmp-tempEnd)*(F_decay)^(f-1)+tempEnd;
    end
    
% ========================================================================
elseif mode == "adaptive"
    disp("'Adaptive' temp_scheduling should be implemented in the solver")
    disp("Setting temp_sched as tempStart (scaler) ...")
    temp_sched = tempStart;
    
% ========================================================================
else
    error("temp_sched_mode doesn't exist")
end



end

