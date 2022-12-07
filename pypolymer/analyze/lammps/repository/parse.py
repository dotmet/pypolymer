import numpy as np
import os
import sys
import matplotlib.pyplot as plt

file_name = 'log.lammps'
num_point = 20
try:
    file_name = sys.argv[1]
except:
    pass

try:
    num_point = sys.argv[2]
except:
    pass
data_matrix = []
start_id = 0
end_id = 0

with open(file_name, 'r') as f:
    raw_data = f.readlines()
    for line in raw_data:
        if "Step" in line:
            start_id = raw_data.index(line, end_id)
        if ("Loop" in line):
            end_id = raw_data.index(line, start_id)
    print(start_id, end_id)
    if start_id >= end_id:
        end_id = len(raw_data)
    data_matrix = raw_data[start_id:end_id+1]
    f.close()

data_matrix =list(map(lambda x:x.split()[:-1], data_matrix))
plot_data = np.array(data_matrix[1:-1], dtype=np.float64)
x_axis = (plot_data[:, 0]-90000)/10000
y_axis = plot_data[:, -1]

#########################
lc='bgrcmykw'
ls=['-','--','-.',':','']
mk='.,ov^<>1234sp*hH+xDd|_'
#########################

idl = int((len(x_axis)-1)/num_point)
new_data = []
for i in range(num_point+1):
    idx = idl*i
    new_data.append([x_axis[idx], y_axis[idx]])
new_data.append([x_axis[-1], y_axis[-1]])
new_data = np.array(new_data)
plt.plot(new_data[:, 0], new_data[:, 1], marker=mk[4], color=lc[2], linestyle=ls[0])
#plt.plot(new_data[:, 0]+100, new_data[:, 1]+1)
plt.legend(["b=300"])
plt.xlabel("Time step(10000" + chr(964)+")")
plt.ylabel("Total Energy")
plt.savefig("300b.png", dpi=500)
#plt.grid()
plt.show()
