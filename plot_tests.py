import timeit
from aks_implementation import aks_prime, fermat_sm_generate_lhs, fermat_bin_generate_lhs
from aks_implementation import rfind, FAST_POWERING_BUFFER
from polynomial import binomial_coefficient
import matplotlib
from math import log

TEST_SIZE = 4000
results = {}
for n in range(3000, TEST_SIZE):
    s = timeit.default_timer()
    for _ in range(5):
        aks_prime(n)
    e = timeit.default_timer()
    results[n] = (e - s)/5

import matplotlib.pyplot as plt

fig = plt.figure()

##exp_y = [pow(1 + 1e-56, _) for _ in range(TEST_SIZE)]
##exp = fig.add_subplot(1,1,1)
##line = exp.plot([_ for _ in range(TEST_SIZE)], exp_y, color='red', lw=1, label = '(1 + 1e-56)^x')
##exp.set_yscale('log')


x, y = list(results.keys()), list(results.values())
output = fig.add_subplot(1,1,1)
line = output.plot(x, y, color='blue', lw=1, label = 'AKS runtime')
output.set_yscale('log')

##log13_y = [log(_, 2)**13 for _ in range(1, TEST_SIZE)]
##log13 = fig.add_subplot(1,1,1)
##line = log13.plot([_ for _ in range(1, TEST_SIZE)], log13_y, color='green', lw=1)
##log13.set_yscale('log')

log1_y = [log(_, 2) for _ in range(1, TEST_SIZE)]
log1 = fig.add_subplot(1,1,1)
line = log1.plot([_ for _ in range(1, TEST_SIZE)], log1_y, color='green', lw=1, label = 'log(x, 2)')
log1.set_yscale('log')

output.set_title('AKS Runtime Comparison', color='black')
output.set_xlabel('Integer', color='black')
output.set_ylabel('Time (seconds)', color='black')

plt.legend(loc = 0)
plt.show()


'''
#saved times, per line:
#(n, time)

TEST_SIZE = 4000
import matplotlib.pyplot as plt
from math import log
f = open('D:/waterloo/co685/aks/code/plot_outputs/draft4/all.txt', 'r')
results = {}
for line in f.readlines():
    n = float(line.strip().split(',')[0][1:])
    time = float(line.strip().split(',')[1][:-1])
    results[n] = time
f.close()
fig = plt.figure()
x, y = list(results.keys()), list(results.values())
output = fig.add_subplot(1,1,1)
line = output.plot(x, y, color='blue', lw=1, label = 'AKS runtime')
output.set_yscale('log')
_110_log12_y = [(1e-10)*log(_, 2)**12 for _ in range(1, TEST_SIZE)]
_110_log12 = fig.add_subplot(1,1,1)
line = _110_log12.plot([_ for _ in range(1, TEST_SIZE)], _110_log12_y, color='green', lw=1, label = '(1e-10){log(x, 2)}^12')
output.set_title('AKS Runtime Comparison', color='black')
output.set_xlabel('Integer', color='black')
output.set_ylabel('Time (seconds)', color='black')
plt.legend(loc = 0)
plt.show()
'''


'''
#plotting
import matplotlib.pyplot as plt
import time
import numpy as np
from scipy.interpolate import spline

    # Create a canvas to place the subgraphs
canvas = plt.figure()
rect = canvas.patch
rect.set_facecolor('white')

# Iterate through the lines and parse them

x_smooth, y_smooth = list(results.keys()), list(results.values())

# Define the matrix of 1x1 to place subplots
# Placing the plot1 on 1x1 matrix, at pos 1
sp1 = canvas.add_subplot(1,1,1, axisbg='w')
sp1.set_yscale('log')
#sp1.plot(x, y, 'red', linewidth=2)
sp1.plot(x_smooth, y_smooth, 'red', linewidth=0.5)

# Colorcode the tick tabs 
sp1.tick_params(axis='x', colors='red')
sp1.tick_params(axis='y', colors='red')

# Colorcode the spine of the graph
sp1.spines['bottom'].set_color('r')
sp1.spines['top'].set_color('r')
sp1.spines['left'].set_color('r')
sp1.spines['right'].set_color('r')

# Put the title and labels
sp1.set_title('AKS Runtimes', color='red')
sp1.set_xlabel('Integer', color='red')
sp1.set_ylabel('Time (seconds)', color='red')

# Show the plot/image
plt.tight_layout()
plt.grid(alpha=0.8)
plt.savefig("example6.eps")
plt.show()
'''
