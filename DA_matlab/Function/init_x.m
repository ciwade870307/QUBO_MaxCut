function [x] = init_x(Problem,Param)

if Param.mode.x_init_mode == "rand"
    x = binornd(1,0.5*ones(Problem.N,1));
elseif Param.mode.x_init_mode == "zeros"
    x = zeros(Problem.N,1);
elseif Param.mode.x_init_mode == "ones"
    x = ones(Problem.N,1);
elseif Param.mode.x_init_mode == "eye"
    X = eye(sqrt(Problem.N));
    x = reshape(X,[Problem.N,1]);
elseif Param.mode.x_init_mode == "spec"
    x = Param.mode.x_init;
else
    error("Incorrect x initialization!")
end
    
end

