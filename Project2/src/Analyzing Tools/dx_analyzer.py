import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import seaborn as sns
from time import time


def read_parameters(folder):
    """
    Take the folder where the parameter file is found, and returns
    a dictionary with the name of the variables as keys,
    and the value of parameter as the value.
    """
    parameters = {}
    with open(folder+"metadata.txt") as f:
        for line in f:
            words = line.split()
            if len(words) == 0:
                continue
            if words[0][0] == "#":
                continue
            parameters[words[0]] = eval(words[1])

    return parameters

def block(x):
    """
    The blocking function made by Marius. It is made so to give an error if
    the number if MC cycles times the number of processes are not on the form
    2^n.
    """
    # preliminaries
    n = len(x);
    if abs(np.log2(n) - int(np.log2(n))) > 1e-5:
        print("Number of MC_cycles are not on the form 2^n. Change that " \
                            +"and come back an other time")

        print(np.log2(n))
        exit()
    d = int(np.log2(n)); s, gamma = np.zeros(d), np.zeros(d);
    mu = np.mean(x); t0 = time()

    # estimate the auto-covariance and variances
    # for each blocking transformation
    for i in np.arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*np.sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = np.var(x)
        # perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2])

    # generate the test observator M_k from the theorem
    M = (np.cumsum( ((gamma/s)**2*2**np.arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # we need a list of magic numbers
    q =np.array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in np.arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")
    ans = s[k]/2**(d-k)
    # print("Runtime: %g sec" % (time()-t0)); print("Blocking Statistics :")
    # print("average            iterations      std. error")
    # print("%8g %20g %15g" % (mu, k, ans**.5))
    return ans


class getVariance:
    """
    Is is an ad hoc made class to get the Evolving variance of different data,
    with different dx(dt). This is in no way standarized, so to get it to work
    the data has to be in files with the same names as in the init.
    """
    def __init__(self,folder):
        self.folder = folder
        self.parameters = read_parameters(folder)
        self.list_of_dx = [1.5,1.25,1,0.75,0.5,0.25,0.1] #For brute force
        #self.list_of_dx = [1,0.5,0.1,0.05,0.01,0.005,0.001] #For importance

        self.file_names = []
        self.accept_file_names = []
        for dx in self.list_of_dx:
            self.file_names.append("data_%1.6f.bin" %(dx))
            self.accept_file_names.append("accept_data_%1.6f.bin" %(dx))
        self.data = np.zeros((self.parameters["MC_cycle"],len(self.list_of_dx)))
        self.accept_data = np.zeros(len(self.list_of_dx))
        self.read_data()
        self.plot_blocking()

    def read_data(self):
        for index,f in enumerate(self.file_names):
            self.data[:,index] = np.fromfile(self.folder + f,sep=" ")
        for index,f in enumerate(self.accept_file_names):
            self.accept_data[index] = np.mean(np.fromfile(self.folder + f,sep=" "))


    def plot_blocking(self):
        """
        Plots the gradient of the energy against dx, and average energy
        against dx. The energy uses blocking to get the error.
        """
        sns.set_style("darkgrid")
        plt.rcParams.update({'font.size':12})
        plt.rcParams['mathtext.fontset']='stix'
        plt.rcParams['font.family']='STIXGeneral'

        fig = plt.figure()
        ax = plt.subplot(111)

        mean_data = np.zeros_like(self.accept_data)
        error_data = np.zeros_like(self.accept_data)

        for i,dx in enumerate(self.list_of_dx):
            mean_data[i] = np.mean(self.data[:,i])

            error_data[i] = block(self.data[:,i])

        plt.plot(self.list_of_dx,self.accept_data,"b.-",markersize=10) #For brute force
        #plt.semilogx(self.list_of_dx,self.accept_data,"b.-",markersize=10) #For importance

        plt.title("Acceptance Rate for Different dx",fontsize=20) #For brute force
        #plt.title("Acceptance Rate for Different dx with Importance Sampling",fontsize=20) #For importance

        plt.xlabel("dx",fontsize=20)
        plt.ylabel("Acceptance Rate",fontsize=20)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlim(self.list_of_dx[0],self.list_of_dx[-1])
        plt.show()


        ax2 = plt.subplot(111)

        ax2.errorbar(self.list_of_dx,mean_data,yerr=error_data,fmt='b.-',markersize=10)

        ax2.get_yaxis().get_major_formatter().set_useOffset(False)
        ax2.get_xaxis().get_major_formatter().set_useOffset(False)
        ax2.set_xlabel("dx",fontsize=20)
        ax2.set_ylabel(r"$\langle E \rangle$",fontsize=20)
        ax2.set_xlim(self.list_of_dx[0],self.list_of_dx[-1])
        plt.title(r"$\langle E \rangle$ for Different dx",fontsize=20) #For brute force
        #plt.title(r"$\langle E \rangle$ for Different dx with Importance Sampling",fontsize=20) #For importance
        plt.tight_layout()
        plt.autoscale()
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.show()






if __name__ == '__main__':
    getVar = getVariance("../output/")
