import time
import numpy as np
import pickle
from glob import glob
from tqdm import tqdm
import sys
sys.path.insert(0, '..')
from DA4 import DA

def bqp_test(path="../data/g05/", file="g05", iters=1000, spin=False, b_matrix=False):
    data_list = [tag for tag in glob("{}/{}*".format(path, file))]
    solution_file = open("{}/solution.txt".format(path), "r").readlines()
    solution = {}
    for s in solution_file:
        f, v = s.split(" ")
        if b_matrix == True :
            solution[f] = int(v)
        else :
            solution[f] = -int(v)

    output_file = open("hist_result/da_{}hist{}.csv".format(file, iters), "w")
    output_file.write("{},{},{},{},{},{},{},{},{}\n".format("Instance","BKS", "Cut_best", "Acc_best", "Cut_avg", "Acc_avg", "Std", "t_avg", "# Opts/ Total runs"))
    data_list = np.sort(data_list)
    run=500
    for data in tqdm(data_list):
        sv = solution[data.split("/")[-1]]
        distribution = np.zeros(10)
        t = 0
        min_e = 0
        f = open(data, "r").readlines()

        bin_size = f[0].split(" ")[0]
        Q = np.zeros([int(bin_size) + 1, int(bin_size) + 1])
        init_bin = np.zeros([int(bin_size) + 1])
        init_bin[-1] = 1
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

        E_Q_list = []
        for n in tqdm(range(run)):
   
            #### annealing ####
            init_e = init_bin @ Q
            da1 = DA(Q, init_bin, maxStep=iters,
                     betaStart=0.01, betaStop=1.61, kernel_dim=(32 * 2,), spin=spin, energy = init_e)
            da1.run()
            # print(da1.binary)
            # print(f'time spent: {da1.time}')

            bin1 = np.expand_dims(da1.binary, axis=1)
            output = np.matmul(np.matmul(bin1.T, Q), bin1)[0][0]
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
            t += da1.time
    
        output_file.write("{},{},{},{:.2f},{:.2f},{:.2f},{:.2f},{:.3f},{}\n".format(data.split("/")[-1].split(".")[0], abs(sv), abs(min_e), abs(100*min_e/sv), 
        abs(np.mean(E_Q_list)), abs(100*np.mean(E_Q_list)/sv), abs(100*np.std(E_Q_list)/sv), t/run, str(distribution[0])+"/"+str(run)))
        output_file.flush()


if __name__ == '__main__':

    # bqp_test(file="g05_100", iters=40000)
    bqp_test(path="../data/gset/", file="G", iters=40000)
    # bqp_test(path="../data/beasley/", file="bqp50-1", iters=40000, b_matrix=True)
