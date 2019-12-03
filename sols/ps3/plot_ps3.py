import matplotlib.pyplot as plt
import numpy


# M = 2
l = [0,1,2,3,4,5,6,7,8]
var1 = [1.3198e+06, 1.0088e+06, 1.1307e+06, 1.1774e+06, 1.2212e+06, 1.2210e+06, 1.2449e+06, 1.2374e+06, 1.2241e+06]
var2 = [1.3198e+06, 3.9112e+05, 1.4528e+04, 3.1918e+03, 1.0249e+03, 3.8172e+02, 2.0106e+02, 7.7466e+01, 2.1480e+01]

del1 = [2.8964e+03, 3.1685e+03, 3.4816e+03, 3.6241e+03, 3.7036e+03, 3.7547e+03, 3.7569e+03, 3.7677e+03, 3.7726e+03]

chk1 = [0.0000e+00,  1.7359e-02, 1.5991e-01, 2.3864e-01, 4.6202e-02, 2.6271e-01, 3.5923e-01, 2.5087e-02, 9.6277e-04]

kur1 = [4.1337e+00, 5.0677e+00, 3.7604e+00, 6.2504e+00, 2.4051e+01, 8.6250e+01, 2.3836e+02, 7.6681e+02, 2.1189e+03]

epss = [1, 2, 5, 10, 20, 50]
Ns = [[5179697, 1985111,   534068,   170576,    47318,    15471,     6906,     3177,     1176,      332,     174,       98,       52],
      [1273420,   488970,  131551,    42773,    11867,     4638,     1960,      710,      195,       72,       52,       19],
      [ 193260,    73735,    19838,     6788,     1875,      513,      147,       59,       28,       19,        7],
      [ 48650,    18554,     4992,     1698,      473,      205,       69,       19,        8],
      [11643,     4523,     1217,      429,      120,       33,       10,        4,        1],
      [1779,     666,      180,       67,       19,        5,        2]]

mlmc = [3.7844e+03, 3.7839e+03 , 3.7913e+03, 3.7713e+03, 3.7713e+03, 3.7272e+03]

"""
# M = 3
l = [0,1,2,3,4,5,6,7]
var1 = [3.8290e+06, 1.0361e+06, 1.1418e+06, 1.2250e+06, 1.1957e+06, 1.2382e+06, 1.2180e+06, 1.2296e+06]
var2 = [3.8290e+06, 1.5255e+06, 1.4275e+04, 1.9826e+03, 3.9989e+02, 1.0887e+02, 4.4435e+01, 1.0228e+01]

del1 = [2.8874e+03, 3.2361e+03, 3.5874e+03, 3.7277e+03, 3.7631e+03, 3.7799e+03, 3.7843e+03, 3.7809e+03]

chk1 = [0.0000e+00, 3.3235e-02, 2.3148e-01, 3.8606e-01, 1.1255e-01, 6.5891e-02, 5.0735e-04, 1.0518e-01]

kur1 = [5.7218e+00, 7.1935e+00, 3.7094e+00, 1.2290e+01, 8.8003e+01, 5.0056e+02, 1.2370e+03, 3.5736e+03]

epss = [1, 2, 5, 10, 20, 50]
Ns = [[16436148, 5989478, 997975, 218057, 37273, 6290, 1356, 318, 100],
      [4038204, 1468972,   244762,    54359,     9299,     1570,      341,  57],
      [660806,   241253,    40198,     8761,     1499,      387,       77,       13],
      [160426,    59861,     9975,     2173,      372,       63,       11],
      [39946,    14645,     2441,      554,       95,       16],
      [6343,     2163,      361,       84,       14,        3]]

mlmc = [3.7837e+03, 3.7861e+03, 3.7827e+03, 3.7745e+03, 3.7755e+03, 3.8010e+03]
"""

"""
plt.figure(figsize=(8, 6))

plt.subplot(3, 2, 1)
plt.plot(l,     numpy.log2(var2),     '*-',  label=r'$P_l$')
plt.plot(l[1:], numpy.log2(var1[1:]), '*--', label=r'$P_l - P_{l-1}$')
plt.xlabel('level $l$')
plt.ylabel(r'$\mathrm{log}_2(\mathrm{variance})$')
plt.legend(loc='lower left', fontsize='x-small')

plt.subplot(3, 2, 2)
plt.plot(l,     numpy.abs(del1),     '*-',  label=r'$P_l$')
plt.xlabel('level $l$')
plt.ylabel(r'Mean')
plt.legend(loc='lower left', fontsize='x-small')

plt.subplot(3, 2, 3)
plt.plot(l[1:], chk1[1:], '*--')
plt.xlabel('level $l$')
plt.ylabel(r'consistency check')
axis = plt.axis(); plt.axis([0, max(l), axis[2], axis[3]])

plt.subplot(3, 2, 4)
plt.plot(l[1:], kur1[1:], '*--')
plt.xlabel('level $l$')
plt.ylabel(r'kurtosis')
axis = plt.axis(); plt.axis([0, max(l), axis[2], axis[3]])

"""

plt.figure(figsize=(12, 6))

styles = ['o--', 'x--', 'd--', '*--', 's--']
plt.subplot(1, 2, 1)
for (eps, N, style) in zip(epss, Ns, styles):
    plt.semilogy(N, style, label=eps)
plt.xlabel('level $l$')
plt.ylabel('$N_l$')
plt.legend(loc='upper right', frameon=True, fontsize='x-small')

eps = numpy.array(epss)
plt.subplot(1, 2, 2)
plt.plot(eps, mlmc, "*--")
plt.xlabel(r'accuracy $\epsilon$')
plt.ylabel(r'MLMC value')
plt.legend(fontsize='x-small')
axis = plt.axis(); plt.axis([min(eps), max(eps), axis[2], axis[3]])


plt.subplots_adjust(wspace=0.3)
plt.show()

    
    
