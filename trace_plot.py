# %%
import numpy as np
import matplotlib.pyplot as plt
import time

time_duration = 0.15

file = 'phase_space.dat'

# Read normal gaussian from results.dat
data = np.genfromtxt(file, delimiter=',')

r = data[:, :1]
print('SHAPE r')
print(np.shape(r))

r_sq = np.power(r, 2) 

d = r_sq[:, 0]
print('SHAPE d')
print(np.shape(d))

c = np.sqrt(d)
print('SHAPE c')
print(np.shape(c))

fig = plt.figure()
ax = fig.add_subplot()
ax.plot( np.linspace(0, 2, num=len(r[:, 0])), r[:, 0], alpha = 0.2, color='blue', marker='.')
    ##time.sleep(time_duration)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
#ax.set_xlim(-0.1, 0.1)
#ax.set_ylim(-0.1, 0.1)
plt.show()


# %%
v = data[:, 1:]
print('SHAPE v')
print(np.shape(v))

v_sq = np.power(v, 2) 

V = v_sq[:, 0]
print('SHAPE V')
print(np.shape(V))

V = np.sqrt(V)
print('SHAPE c')
print(np.shape(V))

fig = plt.figure()
ax = fig.add_subplot()
ax.plot( np.linspace(0, 2, num=len(v[:, 0])), v[:, 0], alpha = 0.2, color='red', marker='.')
    ##time.sleep(time_duration)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
#ax.set_xlim(-0.1, 0.1)
#ax.set_ylim(-0.1, 0.1)
plt.show()
# %%
