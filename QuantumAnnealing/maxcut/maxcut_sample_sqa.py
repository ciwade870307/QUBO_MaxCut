import time
import numpy as np
import pickle
from glob import glob
from tqdm import tqdm
import sys
sys.path.insert(0, '..')

from DA4 import DA
from qubo_gen import *
from SQA import *


def bqp_test(path="../data/g05/", file="g05", iters=1000, spin=False, b_matrix=False):
    M = 4
    T = 0.15

    data_list = [tag for tag in glob("{}/{}*".format(path, file))]
    solution_file = open("{}/solution.txt".format(path), "r").readlines()
    solution = {}
    for s in solution_file:
        f, v = s.split(" ")
        if b_matrix == True :
            solution[f] = int(v)
        else :
            solution[f] = -int(v)

    output_file = open("hist_result/sqa_{}hist{}_M{}T{}.txt".format(file, iters, M, T), "w")
    data_list = np.sort(data_list)
    run=1
    for data in tqdm(data_list):
        sv = solution[data.split("/")[-1]]
        distribution = np.zeros(10)
        t = 0
        min_e = 0
        f = open(data, "r").readlines()

        bin_size = f[0].split(" ")[0]
        Q = np.zeros([int(bin_size) + 0, int(bin_size) + 0])
        init_bin = np.zeros([int(bin_size) + 0])
        # init_bin[-1] = 1
        for ele in f[1:]:
            i, j, v = ele.split()
            if b_matrix == True :
                Q[int(i) - 1, int(j) - 1] += int(v)
                if (int(i) != int(j)):
                    Q[int(j) - 1, int(i) - 1] += int(v)
            else :
                if (int(i) == int(j)) :
                    print('No edge connected at the same Node',int(i),int(j))
                else :
                    Q[int(i) - 1, int(j) - 1] += int(v)
                    Q[int(j) - 1, int(i) - 1] += int(v)
                    Q[int(i) - 1, int(i) - 1] += -int(v)
                    Q[int(j) - 1, int(j) - 1] += -int(v)

        J, h, offset = qubo2ising(Q)
        N = int(bin_size)
        steps = iters
        Gamma0 = 10
        Gamma1 = 1e-8
        decay_rate = (Gamma1 / Gamma0)**(1/(steps-1))
        schedule = [Gamma0 * decay_rate**i for i in range(steps)]
        F = 1
        G = False

        for n in tqdm(range(run)):
            #### annealing ####
            # init_e = init_bin @ Q
            start_time = time.time()
            spin = one_SQA_run(J, h, schedule, M, T, field_cycling=F, return_pauli_z=True, enable_global_move=G)
            total_time = time.time() - start_time
            binary = np.array(spin)>=1
            output = binary.T@Q@binary
            # print(output)
            print("solution", sv, output)

            if output < min_e :
                min_e = output

            if (sv - output) >= 0:
                l = 0
                distribution[0] += 1
                t += total_time
                print('find optimal solution', sv, output)
                break
            else:
                l = int(np.ceil((abs((sv - output) / sv) * 100)))
            if l > 9:
                l = 9
            distribution[l] += 1
            t += total_time
    
        output_file.write("{}\t{}\t{}\t{}\n".format(data.split("/")[-1].split(".")[0], min_e, t, distribution))
        output_file.flush()


if __name__ == '__main__':

    bqp_test(file="g05_100.0", iters=1000)
    # bqp_test(path="../data/gset/", file="G1.txt", iters=1000)
    # bqp_test(path="../data/beasley/", file="bqp50-1.sparse", iters=4000, b_matrix=True)
