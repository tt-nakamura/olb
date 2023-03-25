import numpy as np
import matplotlib.pyplot as plt
from olbers import olbers

a = [0.5, 1, 2]
Omega0 = np.linspace(0.01, 1, 100)
K = 1e-15

plt.figure(figsize=(5, 3.75))

for a in a:
    i = []
    for om in Omega0:
        i.append(olbers(a, om, K))
    plt.semilogy(Omega0, i, label=r'$a=%g$'%a)

plt.axis([0, 1, 1e-17, 1e-13])
plt.legend()
plt.xlabel(r'$\Omega_0$')
plt.ylabel(r'$i = I/I_{\odot}$')
plt.tight_layout()
plt.savefig('fig2.eps')
plt.show()
