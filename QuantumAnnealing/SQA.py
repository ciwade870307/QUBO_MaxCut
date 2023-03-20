from dataclasses import field
import numpy as np
from scipy.sparse import block_diag
import time
import matplotlib.pyplot as plt


def one_SQA_run(J, h, trans_fld_sched, M, T, field_cycling = 1, sd = None, init_state = None, return_pauli_z = False, enable_global_move = False):
    """
    One path-integral Monte Carlo simulated quantum annealing run over the full transverse field strength schedule.
    The goal is to find a state such that sum(J[i, j]*state[i]*state[j]) + sum(h[i]*state[i]) is minimized.
    
    Parameters:
        J (2-D array of float): The matrix representing the coupling field of the problem.
        h (1-D array of float): The vector representing the local field of the problem.
        trans_fld_sched (list[float]): The transeverse field strength schedule for QA.
                                       The number of iterations is implicitly the length of trans_fld_schedule.
        M (int): Number of Trotter replicas. To simulate QA precisely, M should be chosen such that T M / Gamma >> 1.
        T (float): Temperature parameter. Smaller T leads to higher probability of finding ground state.
        sd (default=None): Seed for numpy.random.default_rng().
        init_state (1-D array of int, default=None): The boolean vector representing the initial state.
                                                     If None, a random state is chosen.
        return_pauli_z (bool, default=False): If True, returns a N-spin state averaged over the imaginary time dimension.
                                              If False, returns the raw N*M-spin state.
        enable_global_move (bool, default=Falss): If True, apply global move technique.
        field_cycling (int, default=1): Numbers of cycles in field-cycling. field_cycling=1 is equivalent to not apply field-cycling technique.
    
    Return: final_state (1-D array of int)
    """

    rng = np.random.default_rng(seed = sd)

    # normalize
    norm_coef = np.sqrt(J.shape[0] / (np.sum(J**2) + 0.5 * np.sum(h**2))) # normalization
    J_n = J * norm_coef
    h_n = h * norm_coef


    N = J_n.shape[0]

    J_n = 0.5 * (J_n + J_n.T)
    j = block_diag([J_n / M] * M).toarray()
    h_extended = np.repeat(h_n / M, M)

    Jp_terms = np.eye(N*M, k=N) + np.eye(N*M, k=N*(1-M))
    Jp_terms = 0.5 * (Jp_terms + Jp_terms.T)

    steps = len(trans_fld_sched)
    T_decrement = T / field_cycling

    # Annealing
    for cycle in range(field_cycling):
        if init_state is None:
            state = 2 * rng.binomial(1, 0.5, N*M) - np.ones(N*M)
        else:
            state = np.tile(init_state, M)

        for Gamma in trans_fld_sched[0::field_cycling]:
            Jp_coef = -0.5 * T * np.log(np.tanh(Gamma / (M * T)))

            # Global move
            if enable_global_move:
                candidate = 0
                delta_E_min = -4 * j[0].dot(state) * state[0] - 2 * h_extended[0] * state[0]
                for flip_pos in range(N):
                    delta_E = 0
                    for flip_replica in range(M):
                        flip = N * flip_replica + flip_pos
                        delta_E += -4 * j[flip].dot(state) * state[flip] - 2 * h_extended[flip] * state[flip]
                    if delta_E < delta_E_min:
                        delta_E_min = delta_E
                        candidate = flip_pos

                if rng.binomial(1, np.minimum(np.exp(-delta_E_min/T), 1.)):
                    for flip_replica in range(M):
                        flip = N * flip_replica + candidate
                        state[flip] *= -1

            # Local move
            for flip in range(N*M):
                delta_E = -4 * (j[flip] - Jp_coef * Jp_terms[flip]).dot(state) * state[flip] - 2 * h_extended[flip] * state[flip]
                if rng.binomial(1, np.minimum(np.exp(-delta_E/T), 1.)):
                    state[flip] *= -1

        # Field-cycling
        if cycle != field_cycling - 1:
            T -= T_decrement
            for Gamma in trans_fld_sched[-1::-100 * field_cycling]:
                Jp_coef = -0.5 * T * np.log(np.tanh(Gamma / (M * T)))

                for flip in range(N*M):
                    delta_E = -4 * (j[flip] - Jp_coef * Jp_terms[flip]).dot(state) * state[flip] - 2 * h_extended[flip] * state[flip]
                    if rng.binomial(1, np.minimum(np.exp(-delta_E/T), 1.)):
                        state[flip] *= -1


        state = [int(i) for i in state]
        

    if return_pauli_z:
        final_state = state[:N]
        final_E = np.sum(np.multiply(J,np.outer(final_state,final_state))) + np.dot(h,final_state)
        for i in range(M):
            temp_state = state[i*N:(i+1)*N]
            temp_E = np.sum(np.multiply(J,np.outer(temp_state,temp_state))) + np.dot(h,temp_state)
            if temp_E < final_E:
                final_state = temp_state.copy()
                final_E = temp_E
        return final_state
    else:
        return state


##########################################################################################################


def main():

    sd = 7

    np.random.seed(0)
    num_par = np.random.normal(0, 10, 20)
    N = len(num_par)

    J = np.outer(num_par, num_par)
    np.fill_diagonal(J, 0)
    h = np.zeros(N)

    norm_coef = np.sqrt(J.shape[0] / (np.sum(J**2) + 0.5 * np.sum(h**2))) # normalization
    J = J * norm_coef
    h = h * norm_coef

    M = 20
    T = 0.05

    steps = 1000
    Gamma0 = 10
    Gamma1 = 1e-8
    decay_rate = (Gamma1 / Gamma0)**(1/(steps-1))
    schedule = [Gamma0 * decay_rate**i for i in range(steps)]


    # SQA
    start_time = time.time()
    ans = one_SQA_run(J, h, schedule, M, T, sd=sd, return_pauli_z=True)
    total_time = time.time() - start_time

    E_sqa = np.sum(np.multiply(J,np.outer(ans,ans))) + np.dot(h,ans)

    print("-----simulated quantum annealing-----")
    print(f"final state: {ans}")
    print(f"final energy: {E_sqa}; time: {total_time} s")


    # deterministic solver
    from solver import solver
    start_time_solver = time.time()
    ans_solver = solver(J,h)
    total_time_solver = time.time() - start_time_solver

    E_solver = np.sum(np.multiply(J,np.outer(ans_solver,ans_solver))) + np.dot(h,ans_solver)

    print("-----solver-----")
    print(f"ground state: {ans_solver}")
    print(f"ground energy: {E_solver}; time: {total_time_solver} s")

if __name__ == "__main__":
    main()