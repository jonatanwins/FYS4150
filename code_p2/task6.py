import numpy as np
import matplotlib.pyplot as plt 


# Extract eigeinvectors
n1 = 10
N1 = n1 - 1
n2 = 100
N2 = n2 - 1

v1 = np.zeros((N1 + 2, 3)) # N1 + 2 from adding boundary points (0,0) and (1,0)
v2 = np.zeros((N2 + 2, 3)) # N1 + 2 from adding boundary points (0,0) and (1,0)

j1 = np.zeros((N1 + 2, 3)) # N1 + 2 from adding boundary points (0,0) and (1,0)
j2 = np.zeros((N2 + 2, 3)) # N1 + 2 from adding boundary points (0,0) and (1,0)


#---------- Collecting data from files

with open('data/analytical_eigvec_10.txt', 'r') as file:
    for i,line in enumerate(file):
        word = line.split()
        # Three first columns (eigenvectors) corespont to three lowest eigenvalues
        v1[i+1,:] = np.array([word[0],word[1],word[2]])

with open('data/analytical_eigvec_100.txt', 'r') as file:
    for i,line in enumerate(file):
        word = line.split()
        # Three first columns (eigenvectors) corespont to three lowest eigenvalues
        v2[i+1,:] = np.array([word[0],word[1],word[2]])


with open('data/jacobi_eigvec_10.txt', 'r') as file:
    for i,line in enumerate(file):
        word = line.split()
        # Three first columns (eigenvectors) corespont to three lowest eigenvalues
        j1[i+1,:] = np.array([word[0],word[1],word[2]])

with open('data/jacobi_eigvec_100.txt', 'r') as file:
    for i,line in enumerate(file):
        word = line.split()
        # Three first columns (eigenvectors) corespont to three lowest eigenvalues
        j2[i+1,:] = np.array([word[0],word[1],word[2]])


# first and last element in v1 and v2 is equal to zero -> we have boundary conditions



# -------------------------------- Plotting ------------------------------------
# Define x-axis
x1 = np.linspace(0,1,N1 + 2)
x2 = np.linspace(0,1,N2 + 2)
x1 = np.array([x1,x1,x1]).T # To fit to three eigenvectors 
x2 = np.array([x2,x2,x2]).T # To fit to three eigenvectors 


fig = plt.figure(figsize=(10,7))
fig.suptitle(r'Comparison analytical and numerical eigenvectors $v(\hat{x})$')
ax1 = fig.add_subplot(2,1,1)
ax1.set_title(f'n = {n1}')
ax2 = fig.add_subplot(2,1,2)
ax2.set_title(f'n = {n2}')

ax1.plot(x1,v1, label=['Analytical v_0','Analytical v_1','Analytical v_2'])
ax1.plot(x1,j1, label=['Jacobi v_0','Jacobi v_1','Jacobi v_2'])
ax2.plot(x2,v2, label=['Analytical v_0','Analytical v_1','Analytical v_2'])
ax2.plot(x2,j2, label=['Jacobi v_0','Jacobi v_1','Jacobi v_2'])


#ax1.set_xlabel(r'$\hat{x}$') #more minimalistic plots
ax1.set_ylabel('y')
ax1.legend()

ax2.set_xlabel(r'$\hat{x}$')
ax2.set_ylabel('y')
ax2.legend()
plt.savefig('plot_6.pdf')
plt.tight_layout()
plt.show()


