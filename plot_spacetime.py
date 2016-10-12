#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

filename = "debugout.txt"
spin = 3
spin_colors = ["r", "g", "b", "y", "k"]

# Coding of spin states
if spin == 2:
    n_bits = 1
    mask = 1
elif spin == 3 or 4: 
    n_bits = 2
    mask = 3
elif spin == 5:
    n_bits = 3
    mask = 7
    
def get_op_state(op_state, site):
    return (op_state >> (n_bits*site)) & mask


with open(filename) as f:
    for lineidx, line in enumerate(f):
        ops, spins = line.split("[")
        spins = np.fromstring(spins.replace("]\n", ""), dtype=int, sep=" ")
        print spins
        ops = ops.split(")(")
        ops[0] = ops[0].replace("(", "")
        ops[-1] = ops[-1].replace(")", "")
        if '' in ops:
            ops.remove('')
        ops_arr = []
        for op in ops:
            op = op.split(" ")
            ops_arr.append({"s0": int(op[0]), "s1": int(op[1]), "stp0": int(op[2]), "stp1": int(op[3]),\
                            "state": int(op[4]), "time": float(op[5])})

        n_sites = len(spins)
        if lineidx > -100:
        
            # Plot Spin configuration at t=0
            for s in range(spin):
                print np.where(spins == s)
                positions = np.where(spins == s)
                plt.plot(positions[0], np.zeros_like(positions[0]), "o", markersize=10,
                         color = spin_colors[s], label="Spin "+ str(s) )

            # Plot operators
            for op in ops_arr:
                plt.plot([op["s0"], op["s1"]], [op["time"], op["time"]],
                         linestyle="-", color="black", linewidth=10)

            # Plot worldlines
            for i in range(n_sites):
                plt.plot([i, i], [0, 1], linestyle="--", color="grey")

                # get all operators on the line and order to time
                ops_on_line = [op for op in ops_arr if op["s0"]==i or op["s1"]==i]
                ops_on_line.sort(key=lambda x: x["time"])

                current_state = spins[i]
                prev_time = 0;
                for op in ops_on_line:
                    # Entering left
                    if op["s0"]==i:
                        side = 0
                    elif op["s1"]==i:
                        side = 1
                    print "state ", op["state"], get_op_state(op["state"], 0),\
                        get_op_state(op["state"], 1),\
                        get_op_state(op["state"], 2), get_op_state(op["state"], 3)
                    if get_op_state(op["state"], side) != current_state:
                        print get_op_state(op["state"], side), side, current_state
                    assert get_op_state(op["state"], side) == current_state
                    plt.plot([i, i], [prev_time, op["time"]], color=spin_colors[current_state], linewidth = 5)

                    current_state = get_op_state(op["state"], side + 2)
                    prev_time = op["time"]

                plt.plot([i, i], [prev_time, 1], color=spin_colors[current_state], linewidth = 5)
                print i, ops_on_line 


            plt.xlim([-.5, n_sites + .5])
            plt.ylim([-.05, 1.05])
            plt.legend()
            plt.show()
