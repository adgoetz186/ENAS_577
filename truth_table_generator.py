import numpy as np
import igraph as ig
import copy as cp
import matplotlib.pyplot as plt

# sets font for figure text
plt.rcParams.update({'font.size': 50})

def table_to_latex(header,main,main_cancel):
	entries_in_main = np.size(main[0][0]) // 2
	print("c"*entries_in_main)
	print("\\begin{table}[h!]")
	print("\\centering")
	print("\\begin{tabular}{"+"c"*entries_in_main+ "|"+ "c"*entries_in_main+"}")
	print(header + " \\\\ \\hline")
	for i in range(len(main)):
		string_main = main[i][0].astype(str)
		string_main[cancel_list[i] + entries_in_main] = np.char.add(np.char.add("\cancel{" , string_main[cancel_list[i] + entries_in_main]),"}")
		if i < len(main)-1:
			print(str(string_main)[1:-1].replace("'","").replace(" "," & ").replace("\\\\","\\") + " \\\\")
		else:
			print(str(string_main)[1:-1].replace("'","").replace(" "," & ").replace("\\\\","\\"))
	print("\\end{tabular}")
	print("\\end{table}")

def iterate_state(state,function_list):
	new_state = []
	for fun in range(len(function_list)):
		new_state.append(function_list[fun](state))
	return np.array(new_state)

#Problem 1 a
#A = lambda x:  1-x[1]
#B = lambda x:  x[0]
#function_list = [A,B]
#tf_name_list = ["A","B"]

#Problem 1 b (AND)
A = lambda x:  (1-x[1])*(1-x[0])
B = lambda x:  x[0]
function_list = [A,B]
tf_name_list = ["A","B"]

#Problem 1 b (OR)
#A = lambda x:  min((1-x[1])+(1-x[0]),1)
#B = lambda x:  x[0]
#function_list = [A,B]
#tf_name_list = ["A","B"]

# Problem 2
#A = lambda x:  1 - x[2]
#B = lambda x:  1 - x[0]
#C = lambda x:  1 - x[1]

# Problem 3
# A = c
# B = a
# C = a and ~b
# C = a or ~b
#A = lambda x:  x[2]
#B = lambda x:  x[0]
#C = lambda x:  x[0] * (1-x[1])
#C = lambda x:  min(x[0] + (1-x[1]),1)
#function_list = [A,B,C]
#tf_name_list = ["A","B","C"]
state_list = []
for i in function_list:
	if len(state_list) == 0:
		state_list = [[0],[1]]
	else:
		new_list = []
		for i in range(len(state_list)):
			new_list.append(state_list[i]+[1])
			new_list.append(state_list[i] + [0])
		state_list = new_list
		

transition_dict = {}
for state in state_list:
	transition_dict[str(state)] = iterate_state(state,function_list)
edge_list = []
edge_name_list = []
entry_list_lookup = []
cancel_list = []
truth_table = []
for entry in transition_dict.keys():
	entry_list_lookup.append(str(np.fromstring(entry[1:-1], sep=', ', dtype=int)))
for entry in transition_dict.keys():
	array_entry = np.fromstring(entry[1:-1], sep=', ', dtype=int)
	transition_vals = -1*(array_entry-transition_dict[entry])
	transition_indecies = np.nonzero(transition_vals)
	for i in np.nonzero(transition_vals)[0]:
		new_array = cp.deepcopy(array_entry)
		new_array[i] += transition_vals[i]
		edge_list.append((entry_list_lookup.index(str(array_entry)),entry_list_lookup.index(str(new_array))))
		if transition_vals[i] < 0 and np.size(np.nonzero(transition_vals)) > 1:
			edge_name_list.append(f"$\\tau_{{{tf_name_list[i]}}}$")
		elif transition_vals[i] > 0 and np.size(np.nonzero(transition_vals)) > 1:
			edge_name_list.append(f"$t_{{{tf_name_list[i]}}}$")
		else:
			edge_name_list.append("")
	truth_table.append([np.hstack((transition_dict[entry],array_entry))])
	cancel_list.append(np.nonzero(transition_vals)[0])
	print(truth_table,cancel_list)

visual_style = {}
visual_style["edge_label"] = edge_name_list
visual_style["vertex_size"] = 0.3
visual_style["vertex_label"] = entry_list_lookup
visual_style["edge_font"] = 100
visual_style["edge_color"] = 'grey'
print(edge_name_list)
print("A B C a b c")
header = (str(tf_name_list)[1:-1].replace("'","").replace(","," &") + " & " + str(tf_name_list)[1:-1].replace("'","").replace(","," &").lower())
table_to_latex(header,truth_table,cancel_list)
input()
for i in truth_table:
	print(i[0])
g = ig.Graph(edge_list, directed=True)
fig, ax = plt.subplots()
ig.plot(g, target=ax,**visual_style)
plt.show()
