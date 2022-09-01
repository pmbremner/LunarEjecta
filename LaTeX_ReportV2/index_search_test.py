import sys
import numpy as np

x = np.linspace(0., 10., 1000)


a = float(sys.argv[1])

print(x[a>=x], x[a>=x].size)

i = np.argmax(x[a>=x])

print(i, x[i], x[i+1])

m_imp = np.logspace(-6., 6., 100)
U     = np.linspace(5., 80., 20)

X, Y = np.meshgrid(m_imp, U)
Z = np.zeros(np.shape(X))

print(np.shape(Z[0]))