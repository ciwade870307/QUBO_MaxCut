import time
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as la

inst_list = range(1,21)

output_file = open("optSA_Ghist.csv", "w")
output_file.write("{},{},{},{},{},{},{},{},{}\n".format("Instance","BKS", "Cut_best", "Acc_best", "Cut_avg", "Acc_avg", "Std", "t_avg", "# Opts/ Total runs"))

solution_file = open("../example/solution.txt", "r").readlines()
solution = {}
for s in solution_file:
    f, v = s.split(" ")
    solution[f] = -int(v)

for inst in inst_list:
    filename = "G"+str(inst)+".txt"
    input_file = open(filename, "r")
    sv = solution[filename]
    E_Q_list = []
    flag = 1
    for s in input_file:
        if(s[0]!="#"):
            s_split = s.split()
            if(s_split[0]=='Avg'):
                t = s_split[1]
                t = float(t[4:-1])
            elif(s_split[0]=='avg_en:'):
                pass
            else:
                [E_Q_list.append(int(s_split[0])) for i in range(int(s_split[1]))]
                if (flag):
                    num_opt = int(s_split[1])
                    flag = 0

    run = len(E_Q_list)
    min_e = E_Q_list[0]

    output_file.write("{},{},{},{:.2f},{:.2f},{:.2f},{:.2f},{:.3f},{}\n".format(filename, abs(sv), abs(min_e), abs(100*min_e/sv), 
        abs(np.mean(E_Q_list)), abs(100*np.mean(E_Q_list)/sv), abs(100*np.std(E_Q_list)/sv), t, str(num_opt)+"/"+str(run)))
        