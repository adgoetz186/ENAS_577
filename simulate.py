import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

count = 10000
rvs_test = st.halfnorm.rvs(size = (count,1000))
for col_ind in range(np.shape(rvs_test)[1]-1):
	print(col_ind)
	rvs_test[:,col_ind+1] += rvs_test[:,col_ind]
print(rvs_test)
if np.min(rvs_test,axis=0)[-1] < 1:
	print("error")
else:
	print(rvs_test)
	rvs_test[rvs_test <= 1] = 1
	rvs_test[rvs_test>1] = 0
	print(rvs_test)
	
	print(rvs_test)
	print(np.sum(rvs_test)/count)