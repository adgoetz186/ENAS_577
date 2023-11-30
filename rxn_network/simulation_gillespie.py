import rxn_model
import matplotlib.pyplot as plt
import numpy as np
import json
import pandas as pd
import seaborn as sns

# sets font for figure text
plt.rcParams.update({'font.size': 15})

# obtained from https://stackoverflow.com/questions/12309269/how-do-i-write-json-data-to-a-file
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

p_error_B_mean = np.zeros((10,10))
p_error_A_mean = np.zeros((10,10))

p_error_B_var = np.zeros((10,10))
p_error_A_var = np.zeros((10,10))

i_ind = 0
for i in np.logspace(-2,2,10):
    j_ind = 0
    for j in np.logspace(-2,2,10):
        reaction_model = rxn_model.rxn_model()
        reaction_model.add_rxn("bind_A",{"DNA_A_e":1,"A":1},{"DNA_A_f":1},"A*DNA_A_e*kA_on")
        reaction_model.add_rxn("unbind_A",{"DNA_A_f":1},{"DNA_A_e":1,"A":1},"DNA_A_f*kA_off")
        
        reaction_model.add_rxn("bind_B",{"DNA_B_e":1,"B":1},{"DNA_B_f":1},"B*DNA_B_e*kB_on")
        reaction_model.add_rxn("unbind_B",{"DNA_B_f":1},{"DNA_B_e":1,"B":1},"DNA_B_f*kB_off")
        
        reaction_model.add_param({"kA_on":i,"kA_off":1,"kB_on":j,"kB_off":1})
        n = 10000
        tvals = np.linspace(0,150,150)
        result_dict = reaction_model.Gillespie_simulate(tvals,ivc={"DNA_A_e": 1,"DNA_B_e": 1,"A":1,"B":1},n = n)

        #plt.plot(tvals,result_dict["DNA_A_f"][0],label = "DNA_A bound" , color = 'b')
        #plt.plot(tvals, result_dict["DNA_B_f"][0],label = "DNA_B bound", color = 'r')
        #plt.xlabel("time")
        #plt.title(f"q_A = {np.round(i,3)}\nq_B = {np.round(j,3)}")
        #plt.legend()
        #plt.ylabel("abundance")
        #plt.show()
        #plt.plot(tvals, np.average(result_dict["DNA_A_f"],axis = 0), label="DNA_A bound (avg, n = 100)", color='b',alpha = 0.5)
        #plt.plot(tvals, np.average(result_dict["DNA_B_f"],axis = 0), label="DNA_B bound (avg, n = 100)", color='r',alpha = 0.5)
        #plt.title(f"q_A = {np.round(i,3)}\nq_B = {np.round(j,3)}")
        #plt.xlabel("time")
        #plt.legend()
        #plt.ylabel("abundance")
        #plt.show()

        p_error_B_mean[i_ind,j_ind] = np.abs(np.average(result_dict["DNA_B_f"][:,-1]) - j/(j+1) / (j/(j+1)))
        p_error_A_mean[i_ind,j_ind] = np.abs(np.average(result_dict["DNA_A_f"][:,-1]) - i/(i+1) / (i/(i+1)))

        p_error_B_var[i_ind,j_ind] = np.abs((j / (j + 1)**2 - np.var(result_dict["DNA_B_f"][:, -1]))/((j / (j + 1)**2)))
        p_error_A_var[i_ind,j_ind] = np.abs((i / (i + 1)**2 - np.var(result_dict["DNA_A_f"][:, -1]))/((i / (i + 1)**2)))
        print([i,j])
        j_ind +=1
    i_ind += 1
sns.heatmap(p_error_B_mean,vmin = 0,vmax = 1)
plt.title(f"Error in approximating B binding mean\nn={n}")
plt.ylabel("$q_A$")
plt.xticks(np.linspace(1,9,5),np.round(np.logspace(-2,2,5),2))
plt.yticks(np.linspace(1,9,5),np.round(np.logspace(-2,2,5),2))
plt.xlabel("$q_B$")
plt.tight_layout()
plt.show()
sns.heatmap(p_error_A_mean,vmin = 0,vmax = 1)
plt.title(f"Error in approximating A binding mean\nn={n}")
plt.xticks(np.linspace(1,9,5),np.round(np.logspace(-2,2,5),2))
plt.yticks(np.linspace(1,9,5),np.round(np.logspace(-2,2,5),2))
plt.ylabel("$q_A$")
plt.xlabel("$q_B$")
plt.tight_layout()
plt.show()
sns.heatmap(p_error_B_var,vmin = 0,vmax = 1)
plt.title(f"Error in approximating B binding variance\nn={n}")
plt.xticks(np.linspace(1,9,5),np.round(np.logspace(-2,2,5),2))
plt.yticks(np.linspace(1,9,5),np.round(np.logspace(-2,2,5),2))
plt.ylabel("$q_A$")
plt.xlabel("$q_B$")
plt.tight_layout()
plt.show()
sns.heatmap(p_error_A_var,vmin = 0,vmax = 1)
plt.title(f"Error in approximating A binding variance\nn={n}")
plt.xticks(np.linspace(1,9,5),np.round(np.logspace(-2,2,5),2))
plt.yticks(np.linspace(1,9,5),np.round(np.logspace(-2,2,5),2))
plt.ylabel("$q_A$")
plt.xlabel("$q_B$")
plt.tight_layout()
plt.show()