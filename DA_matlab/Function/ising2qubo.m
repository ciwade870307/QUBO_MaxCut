function [Q, Q_offset] = ising2qubo(J,h)
% Convert a Ising problem into a QUBO problem matrix

J = (J+J')/2; % making sure Q is symmetric, not necessary in this file
Q = 4*J + diag(2*h - 4*sum(J,1)');
Q_offset = sum(J,'all') - sum(h);

end