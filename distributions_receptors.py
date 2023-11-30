import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

# sets font for figure text
plt.rcParams.update({'font.size': 15})

x = np.arange(20)
poisson_vals = st.poisson.pmf(x,4)
binom_vals = st.binom.pmf(x,8,.5)

norm_vals_cont = st.norm.pdf(x,4,scale=np.sqrt(2))
norm_vals_low = st.norm.cdf(x-0.5,4,scale=np.sqrt(2))
norm_vals_high = st.norm.cdf(x+0.5,4,scale=np.sqrt(2))
norm_vals_des = norm_vals_high-norm_vals_low
print(norm_vals_des)
print(binom_vals)
print(poisson_vals)

#print(np.sum((x-4)**2*poisson_vals))
#print(np.sum((x-4)**2*binom_vals))
#print(norm_vals)
#print(poisson_vals)
plt.scatter(x,binom_vals,label = "binom",color = 'r')
plt.scatter(x,poisson_vals,label = "poisson", color = 'g')
plt.scatter(x,norm_vals_des,label = "norm (discrete)",color='b')
plt.plot(x,norm_vals_cont,label = "norm (cont)",color='b')
plt.ylabel("pmf/pdf")
plt.xlabel("occupied receptor count")
plt.legend()
plt.tight_layout()
plt.show()