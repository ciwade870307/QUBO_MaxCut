function Problem = run_Problem_norm(Problem)

Problem.Q = Problem.C;
Q_norm_coef = (sqrt(norm(Problem.Q,'fro')));
% Q_norm_coef = 1;
Problem.Q_norm_coef = Q_norm_coef;
Problem.Q = Problem.Q/Q_norm_coef;

end