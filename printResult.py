import numpy as np
#a=np.random.random((5,5))
#np.savetxt('test.txt',a,fmt='%.18e',delimiter=' ')

fname_phi1='/Users/changing/Documents/MATLAB/Graduation/time_50_phi.txt'
fname_phi2='/Users/changing/Documents/MATLAB/Graduation/time_1000_phi.txt'
fname_phi3='/Users/changing/Documents/MATLAB/Graduation/time_2000_phi.txt'
fname_phi4='/Users/changing/Documents/MATLAB/Graduation/time_4000_phi.txt'


a=np.loadtxt(fname_phi1,dtype=np.float64)
a=a.reshape((300,300))

b=np.loadtxt(fname_phi2,dtype=np.float64)
b=b.reshape((300,300))

c=np.loadtxt(fname_phi3,dtype=np.float64)
c=c.reshape((300,300))

d=np.loadtxt(fname_phi4,dtype=np.float64)
d=d.reshape((300,300))


import matplotlib.pyplot as plt
plt.figure()
plt.subplot(221)
plt.imshow(a)
plt.subplot(222)
plt.imshow(b)
plt.subplot(223)
plt.imshow(c)
plt.subplot(224)
plt.imshow(d)

plt.show()
