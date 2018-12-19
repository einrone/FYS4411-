import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

"""
Reads the gradient files and plots the gradient squared of the energy against
iterations for different learning rates and hidden nodes.
"""



sns.set_style("darkgrid")
plt.rcParams.update({'font.size':12})
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['legend.fontsize'] = 16



rates = [0.320000,0.310000,0.300000,0.200000,0.100000]
Ns = [1,2,3,4,5]

colors = ["r","g","b","y", "c"]
ticks = ["*","s",".","8","D"]

fig1=plt.figure()
ax1=fig1.add_subplot(111)
fig2=plt.figure()
ax2=fig2.add_subplot(111)

for rate in range(len(rates)):
    for N in range(len(Ns)):
        filename = "../output/gradient_data_%s_%.6f" %(Ns[N],rates[rate])
        grad = np.fromfile(filename,sep=" ")
        filename_energy = "../output/energy_data_%s_%.6f" %(Ns[N],rates[rate])
        energy = np.fromfile(filename_energy,sep=" ")
        iterations = np.arange(1,len(grad)+1)
        ax1.plot(iterations,np.log(grad),color=colors[rate],marker=ticks[N],label="N=%s;Rate=%s" %(Ns[N],rates[rate]))
        ax2.plot(iterations,energy,color=colors[rate],marker=ticks[N],label="N=%s;Rate=%s" %(Ns[N],rates[rate]))

ax1.legend(prop={'size':30})
ax2.legend(prop={'size':30})
ax1.set_title("Convergence for Different Ns and Learning Rates",fontsize=40)
ax1.set_xlabel("Iterations",fontsize=40)
ax2.set_title("Convergence for Different Ns and Learning Rates",fontsize=40)
ax2.set_xlabel("Iterations",fontsize=40)
ax1.set_ylabel(r"$\log\ |\nabla E_L|^2$",fontsize=40)
ax2.set_ylabel(r"$\langle E \rangle$",fontsize=40)
#ax1.set_ylim(0,10)

box = ax1.get_position()
ax1.set_position([box.x0, box.y0, box.width * 0.92, box.height])

# Put a legend to the right of the current axis
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

box = ax2.get_position()
ax2.set_position([box.x0, box.y0, box.width * 0.92, box.height])

# Put a legend to the right of the current axis
ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
