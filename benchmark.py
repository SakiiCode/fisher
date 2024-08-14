import fisher
import time
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
stats = importr('stats')
import numpy as np


tests = [
	[[8, 8, 3, 5, 2], [5, 3, 3, 0, 2], [8, 9, 9, 0, 0], [9, 4, 5, 3, 2], [4, 6, 6, 1, 0]],
	[[24, 21, 12, 6, 7], [3, 2, 2, 1, 1], [11, 11, 6, 2, 2], [9, 7, 4, 2, 4], [8, 7, 8, 2, 2]]
]

print("-- EXACT TEST --")
for a in tests:
	m = np.array(a)

	start = time.time()
	result = fisher.exact(a,200000000)
	end = time.time()
	print("fisher-rxc", "{:.4f}".format(result), "in {:.2f}s".format(end-start),sep="\t");


	start = time.time()
	result = stats.fisher_test(m, workspace=2e8)[0][0]
	end = time.time()
	print("rpy2\t", "{:.4f}".format(result), "in {:.2f}s".format(end-start),sep="\t")

	print("-----------")
  
print("-- MONTE-CARLO SIMULATION --")
for a in tests:
	m = np.array(a)

	start = time.time()
	result = fisher.sim(a,10000000)
	end = time.time()
	print("fisher-rxc", "{:.4f}".format(result), "in {:.2f}s".format(end-start), sep="\t");


	start = time.time()
	result = stats.fisher_test(m, simulate_p_value=True, B=10000000)[0][0]
	end = time.time()
	print("rpy2\t", "{:.4f}".format(result), "in {:.2f}s".format(end-start), sep="\t")

	print("-----------")