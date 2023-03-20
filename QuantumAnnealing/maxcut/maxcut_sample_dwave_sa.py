#!/usr/bin/env python
# coding: utf-8

# In[1]:


import time
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as la
import pickle
from glob import glob
from tqdm import tqdm
import sys
import neal
import dimod
import math
sys.path.insert(0, '..')
from DA4 import DA
# from qap2qubo import *
# from penalty_weight import *

def bqp_test(path="../data/g05/", file="g05", iters=1000, spin=False, betaStart=0.01, betaStop=0.1, E_offset_increase_rate=0, b_matrix=False):

    sampler = neal.SimulatedAnnealingSampler()
    # print(sampler.parameters['beta_schedule_type'])
    # h = {'a': 0.0, 'b': 0.0, 'c': 0.0}
    # J = {('a', 'b'): 1.0, ('b', 'c'): 1.0, ('a', 'c'): 1.0}
    # sampleset = sampler.sample_ising(h, J, num_reads=10)
    # print(sampleset.first.energy)

    data_list = [tag for tag in glob("{}/{}*".format(path, file))]
    solution_file = open("{}/solution.txt".format(path), "r").readlines()
    solution = {}
    for s in solution_file:
        f, v = s.split(" ")
        if b_matrix == True :
            solution[f] = int(v)
        else :
            solution[f] = -int(v)

    output_file = open("hist_result/dwaveSA_{}hist{}.csv".format(file, iters), "w")
    # output_file = open("hist_result/test", "w")
    output_file.write("{},{},{},{},{},{},{},{},{}\n".format("Instance","BKS", "Cut_best", "Acc_best", "Cut_avg", "Acc_avg", "Std", "t_avg", "# Opts/ Total runs"))
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
        init_bin = np.zeros((int(bin_size),1), dtype=np.bool)
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

        Q_dict = {}
        linear_dict = {}
        quadratic_dict = {}
        for i in range(int(bin_size)):
            for j in range(int(bin_size)):
                # Q_dict['('+str(i+1)+','+str(j+1)+')'] = Q[i,j]
                Q_dict[i,j] = Q[i,j]
                if i==j:
                    linear_dict[i] = Q[i,j]
                else:
                    quadratic_dict[i,j] = Q[i,j]

        # print(linear_dict)
        bqm = dimod.BinaryQuadraticModel(linear_dict, quadratic_dict, 0, dimod.BINARY)
        E_Q_list = []
        for n in tqdm(range(run)):
            #### annealing ####
            # init_e = init_bin @ Q
            start_time = time.time()
            # x = myDA(Q, init_bin, maxStep=iters, betaStart=betaStart, betaStop=betaStop, E_offset_increase_rate=E_offset_increase_rate)
            # sampleset = sampler.sample_qubo(Q_dict,
            #                                 # beta_range=[0.01, 1.61],
            #                                 num_sweeps=iters
            #                                 # beta_schedule_type='geometric'
            #                                 )
            sampleset = sampler.sample(bqm,
                                        # beta_range=[0.01, 1.61],
                                        num_sweeps=iters
                                        # beta_schedule_type='geometric'
                                        )
            # print(sampleset.first.energy)
            # [hot_beta, cold_beta] = neal.default_beta_range(bqm)
            # print(f"hot_beta:{hot_beta}, cold_beta:{cold_beta}")
            total_time = time.time() - start_time
            # output = np.squeeze(x.T@Q@x)
            output = sampleset.first.energy
            # print(output)
            print("solution", sv, output)
            E_Q_list.append(output)


            if output < min_e :
                min_e = output

            if (sv - output) >= 0:
                l = 0
                distribution[0] += 1
                # t += total_time
                print('find optimal solution', sv, output)
                # break
            else:
                l = int(np.ceil((abs((sv - output) / sv) * 100)))
                if l > 9:
                    l = 9
                distribution[l] += 1
            t += total_time
    
        output_file.write("{},{},{},{:.2f},{:.2f},{:.2f},{:.2f},{:.3f},{}\n".format(data.split("/")[-1].split(".")[0], abs(sv), abs(min_e), abs(100*min_e/sv), 
        abs(np.mean(E_Q_list)), abs(100*np.mean(E_Q_list)/sv), abs(100*np.std(E_Q_list)/sv), t/run, str(distribution[0])+"/"+str(run)))
        # output_file.write("{}\t{}\t{}\t{}\n".format(data.split("/")[-1].split(".")[0], min_e, t, distribution))
        output_file.flush()

# bqp_test(file="g05_100", iters=40000, betaStart=0.1, betaStop=10000, E_offset_increase_rate=0.1)
bqp_test(path="../data/gset/", file="G", iters=1000)
# bqp_test(path="../data/gset/", file="G.txt", iters=400)