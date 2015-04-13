#!/usr/bin/env python

import numpy as np
from scipy.signal import lfilter, tf2zpk
import pdb
import matplotlib.pyplot as plt
import pdb

# Sample transfer function H:
n  = 400;
f  = np.logspace(np.log10(1), np.log10(1000), n);  # sample with logarithmic frequency scale
w  = 2*np.pi*f;
R  = 3000
L  = 100
C  = 0.098e-6
Rd = 30000;

H = 1 + 1/Rd* 1. / ( (R + 1j*w*L)**-1 + 1j*w*C);


print H
pdb.set_trace()

b = np.array([1.0770389918072, 2.03115724017531, 0.961079950480885])
a = np.array([1.0000000000000, 1.98606828810359, 0.992414110487639])


z, p, k = tf2zpk(b, a)

plt.figure()
sp = plt.subplot(1,1,1, aspect='equal')
sp.tick_params(axis='both',  labelsize=7)
plt.title("Pole-Zero plot", fontsize=8)
plt.xlabel('Re', fontsize=4)
plt.ylabel('Im', fontsize=4)

plt.plot(np.real(z), np.imag(z), '.k', markersize=8)
plt.plot(np.real(p), np.imag(p), 'xk', markersize=8)

circ = plt.Circle((0,0), radius=1, color='g', fill=False)
sp.add_patch(circ)
plt.xlim([-1.5, 1.5])
plt.ylim([-1.5, 1.5])
#plt.grid()


# Impulse Response of filter
inp  = np.zeros(firLen)

pdb.set_trace()

shn  = lfilter(b,a, inp)


#plt.show()

#pdb.set_trace()


