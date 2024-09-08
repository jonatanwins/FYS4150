import numpy as np 
import matplotlib.pyplot as plt 

# -----------------------Task 2-----------------------------
# Extract data from data/output.txt
x = []
u = []
with open('data/output.txt', 'r') as file:
    for line in file:
        ob = line.split()
        x.append(float(ob[0]))
        u.append(float(ob[1]))
        


fig, ax = plt.subplots()
ax.plot(x,u,label='u(x)')
ax.set_xlabel('x')
ax.set_ylabel('u(x)')
ax.set_title(r'$u(x) = 1 - (1 - e^{-10})x - e^{-10x}$')
ax.legend()
plt.savefig('plots/plot_t2.pdf')
plt.show()



#----------------------------Extract data task 7-8----------------------
n = []
x_list = []
u_list = []
v_list = []
data_list = ['data/v_and_u_n_10.txt', 'data/v_and_u_n_100.txt', 'data/v_and_u_n_1000.txt', 'data/v_and_u_n_10000.txt','data/v_and_u_n_100000.txt', 'data/v_and_u_n_1000000.txt', 'data/v_and_u_n_10000000.txt'] 
#Extract data from data/output.txt
for data in data_list:
    with open(data, 'r') as file:
        line_count = sum(1 for line in file)
        # subtract 3 from: len(v*) = len(v) + 2 + 1
        # given endpoints and steps = n + 1
        n.append(line_count - 3) 
        file.seek(0)  # Rewind file pointer back to the beginning
        x = np.linspace(0,1,line_count) #np.zeros(line_count)
        u = np.zeros(line_count)
        v = np.zeros(line_count)
        for i, line in enumerate(file):
            ob = line.split()
            #x[i] = (float(ob[0]))
            v[i] = (float(ob[0]))
            u[i] = (float(ob[1]))   
        x_list.append(x)
        u_list.append(u)
        v_list.append(v)





#-----------------------------Task7--------------------------------------
fig, ax = plt.subplots()
fig.suptitle('Comparison exact solution, u(x), against numeric solution, v(x)')
ax.plot(x_list[-1],u_list[-1], label=r'$u(x) = 1 - (1 - e^{-10})x - e^{-10x}$') #plotting analytical solution
for i in range(len(n)):
    ax.plot(x_list[i],v_list[i], label=fr'v(x), n_{{steps}} = {n[i]}') #plotting numerical solutions

ax.set_xlabel(r'x_i')
ax.set_ylabel('y')
ax.legend()
plt.savefig('plots/plot_7.pdf')
plt.show()




#-----------------------------Task8--------------------------------------


max_logRel = [] # Task 8c, plot comparison max_relativ error vs n_steps
fig = plt.figure(figsize=(10,10))
fig.suptitle('Errors', fontsize=16)
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
for i in range(len(n)):
    # Not include first and last element (v_ and u_ = 0)
    u_ = u_list[i][1:-2]
    v_ = v_list[i][1:-2]
    x_ = x_list[i][1:-2]

    # ------------- Log - Abs. Errror  -----------------
    logDelta = np.log10(np.abs(u_-v_))
    #print(logDelta)

    ax1.plot(x_,logDelta, label=f'n = {n[i]}')

    # ------------- Log - Relativ Error -----------------
    logRel = np.log10(np.abs( (u_-v_) / u_ ))
    #print(logRel)
    max_logRel.append(max(logRel))

    ax2.plot(x_,logRel, label=f'n = {n[i]}')

ax1.set_title(r'Logarithm of the absolute error,  $log_{10}(|u_i - v_i|) = log_{10}(\Delta_i)$')
ax1.set_xlabel(r'$x_i$')
ax1.set_ylabel(r'$log_{10}(\Delta_i)$')
ax1.legend()
ax1.grid()

ax2.set_title(r'Logarithm of the relative error,  $log_{10}(|\frac{u_i - v_i}{u_i}|) = log_{10}(\epsilon_i)$')
ax2.set_xlabel(r'$x_i$')
ax2.set_ylabel(r'$log_{10}(\epsilon_i)$')
ax2.legend()
ax2.grid()
plt.savefig('plots/plot_8.pdf')
plt.show()


fig = plt.figure(figsize=(8,8))
fig.suptitle(r'$Maximum Relativ Error VS log_{10}(n_{steps})$')
ax = fig.add_subplot(1,1,1)
ax.plot(np.log10(n),max_logRel)
ax.plot(np.log10(n),max_logRel, 'bo')
ax.set_xlabel(r'$log_{10}(n_{steps})$')
ax.set_ylabel(r'$Max(\epsilon_i)$')
ax.grid()
plt.savefig('plots/plot_8_max.pdf')

# -----------------------Task 10-----------------------------

n = []
runtime = []
runtime_special = []
with open('data/runtime', 'r') as file:
    for line in file:
        ob = line.split()
        n.append(float(ob[0]))
        runtime.append(float(ob[1]))
        runtime_special.append(float(ob[2]))
        
fig, ax = plt.subplots()
ax.plot(n, runtime, linestyle='-', marker='*', label='General algorithm')
ax.plot(n, runtime_special, linestyle='-', marker='o', label='Special algorithm')
ax.set_xscale('log')
ax.set_xlabel(r'$\log{n_{steps}}$')
ax.set_ylabel('Runtime [s]')
ax.set_title('Runtime')
plt.savefig('plots/plot_t10.pdf')
plt.legend()
plt.show()
