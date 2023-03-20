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
sys.path.insert(0, '..')
from DA4 import DA
# from qap2qubo import *
# from penalty_weight import *


# In[2]:


def myDA(Q, init_bin, maxStep=1000, betaStart=0.1, betaStop=50, E_offset_increase_rate=0, plot_p=True):
    
    N = Q.shape[0]
    Q_norm_coef = la.norm(Q, 'fro')*N
    Q = Q/np.sqrt(Q_norm_coef)

    x = init_bin
    init_e = Q@init_bin
    
    # betaList = np.linspace(betaStart, betaStop, maxStep)
    tempStart = 1/betaStart
    tempEnd = 1/betaStop
    decay = (tempEnd / tempStart)**(1/(maxStep-1))
    betaList = [1/(tempStart * decay**i) for i in range(maxStep)]

    E_offset = 0
    E_best = 999999999
    pList = []

    for idx_step, beta in enumerate(betaList):
        # I: Trial phase
        idx_ones = np.squeeze(x)
        sum_partial = np.sum(Q[:,idx_ones],1).reshape(-1,1)
        dE = np.multiply((2*(1-2*x)),sum_partial)+np.diag(Q).reshape(-1,1)
        # p = np.exp(-(dE-E_offset)*beta)
        # idx_rand = np.argmax(p)
        # print(idx_rand)
        p = np.minimum(np.exp(-(dE-E_offset)*beta), 1.) # Acceptance Probility 
        accepted = np.random.binomial(1, np.squeeze(p), N)
        # II: Update phase
        if( np.any(accepted) ):
        # if( ~idx_rand ):
            idx_rand = np.random.choice(accepted.nonzero()[0])
            x[idx_rand] ^= True
            E_offset = 0
            pList.append(p[idx_rand])
        else:
            E_offset += E_offset_increase_rate
            pList.append(0)

        E_current = x.T@Q@x
        if( E_current < E_best):
            E_best = E_current
            x_best = x 

    if plot_p:
        plt.figure(figsize=(8,4))
        plt.grid()
        plt.plot(pList, label="Acceptance Prob.", color='r', marker='.', linewidth=1.0, linestyle='')
        # plt.plot(1/np.array(betaList), label="Acceptance Prob.", color='r', marker='', linewidth=1.0, linestyle='-')
        plt.legend(loc='upper right')
        # plt.xticks(self.N*self.self.maxStep)
        plt.xlabel('MC steps', fontsize=13)
        plt.ylabel('Acceptance Prob.', fontsize=13)
        plt.show()

    return x_best


# In[3]:


def bqp_test(path="../data/g05/", file="g05", iters=1000, spin=False, betaStart=0.01, betaStop=0.1, E_offset_increase_rate=0, b_matrix=False):

    data_list = [tag for tag in glob("{}/{}*".format(path, file))]
    solution_file = open("{}/solution.txt".format(path), "r").readlines()
    solution = {}
    for s in solution_file:
        f, v = s.split(" ")
        if b_matrix == True :
            solution[f] = int(v)
        else :
            solution[f] = -int(v)

    output_file = open("hist_result/myda_{}hist{}.txt".format(file, iters), "w")
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

        for n in tqdm(range(run)):
            #### annealing ####
            # init_e = init_bin @ Q
            start_time = time.time()
            x = myDA(Q, init_bin, maxStep=iters, betaStart=betaStart, betaStop=betaStop, E_offset_increase_rate=E_offset_increase_rate)
            total_time = time.time() - start_time
            output = np.squeeze(x.T@Q@x)
            # print(output)
            print("solution", sv, output)

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
    
        output_file.write("{}\t{}\t{}\t{}\n".format(data.split("/")[-1].split(".")[0], min_e, t, distribution))
        output_file.flush()


# In[4]:


bqp_test(file="g05_100.0", iters=10000, betaStart=0.1, betaStop=10000, E_offset_increase_rate=0.1)

