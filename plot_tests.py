import timeit
from aks_implementation import aks_prime
import matplotlib
from math import log

TEST_SIZE = 500
results = {}
for n in range(1, TEST_SIZE):
    s = timeit.default_timer()
    for _ in range(5):
        aks_prime(n)
    e = timeit.default_timer()
    results[n] = (e - s)/5

import matplotlib.pyplot as plt

fig = plt.figure()

x, y = list(results.keys()), list(results.values())
output = fig.add_subplot(1,1,1)
line = output.plot(x, y, color='blue', lw=1, label = 'AKS runtime')
output.set_yscale('log')

log1_y = [log(_, 2) for _ in range(1, TEST_SIZE)]
log1 = fig.add_subplot(1,1,1)
line = log1.plot([_ for _ in range(1, TEST_SIZE)], log1_y, color='green', lw=1, label = 'log(x, 2)')
log1.set_yscale('log')

output.set_title('AKS Runtime Comparison', color='black')
output.set_xlabel('Integer', color='black')
output.set_ylabel('Time (seconds)', color='black')

plt.legend(loc = 0)
plt.show()
