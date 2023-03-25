import numpy as np
import matplotlib.pyplot as plt
from olbers import CosmicTime,olbers

K = 1e-15
Omega0 = [1, 1, 0.2, 0.2]
q1 = [0, 2, 0, 2]
q2 = [np.inf, 20, np.inf, 20]
color = ['r', 'r--', 'b', 'b--']
a = np.linspace(1e-3, 2, 200)

plt.figure(figsize=(5, 3.75))

for om,q1,q2,c in zip(Omega0, q1, q2, color):
    t = CosmicTime(a, om)
    t0 = CosmicTime(1, om)
    i = olbers(a, om, K, q1, q2)
    i0 = olbers(1, om, K, q1, q2)
    l = r'$(\Omega_0,q_1,q_2)=(%g,%g,%g)$'%(om,q1,q2)
    plt.semilogy(t, i, c, label=l.replace('inf', '\infty'))
    plt.semilogy([t0], [i0], 'k*')

i = -np.expm1(-K*t)# static universe
plt.semilogy(t, i, 'k:', label='static')

plt.axis([0, 1.5, 1e-17, 1e-13])
plt.legend()
plt.xlabel(r'$H_0t$')
plt.ylabel(r'$i = I/I_{\odot}$')
plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()
