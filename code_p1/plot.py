import numpy as np 
import matplotlib.pyplot as plt 


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
plt.savefig('plot_t1.pdf')
plt.show()