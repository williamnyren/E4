# %%
import numpy as np
import matplotlib.pyplot as plt
import time

time_duration = 0.15

file = 'phase_space.dat'
with open(file) as f:
    lines = f.read() ##Assume the sample file has 3 lines
    str_settings = lines.split('\n', 1)[0]
f.close()
setting = str_settings.split(",")
dt = float(setting[0])
n_timesteps = int(setting[1])

print(dt, n_timesteps)
# Read normal gaussian from results.dat
data = np.genfromtxt(file, delimiter=',', skip_header=1)


r = data[:n_timesteps, :1]

fig = plt.figure()
ax = fig.add_subplot()
ax.plot( np.linspace(0, 2, num=len(r[:, 0])), r[:, 0], 
alpha = 0.2, color='blue', marker='.')
##time.sleep(time_duration)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
#ax.set_xlim(-0.1, 0.1)
#ax.set_ylim(-0.1, 0.1)
plt.show()


# %%
print('SHAPE data')
print(np.shape(data))
v = data[:n_timesteps, 1:]
print('SHAPE v')
print(np.shape(v))

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
time_step = 50*0.001
ps = np.abs(np.fft.fft(v))**2

freqs = np.fft.fftfreq(v.size, time_step)
idx = np.argsort(freqs)

plt.plot(freqs[idx], ps[idx])


# %%

v_samples = np.array_split(data[:, 1:], n_timesteps)

time_step = 50*0.001
ps0 = np.abs(np.fft.fft(v_samples[0]))**2
freqs0 = np.fft.fftfreq(v_samples[0].size, time_step)
idx = np.argsort(freqs0)
freqs0 = freqs0[idx].T
ps0 = ps0[idx].T


# %%
for v in v_samples[1:]:
    ps = np.abs(np.fft.fft(v))**2
    freqs = np.fft.fftfreq(v.size, time_step)
    idx = np.argsort(freqs)
    ps0 = np.concatenate((ps0, ps[idx].T), axis=None)
    freqs0 = np.concatenate((freqs0, freqs[idx].T), axis=None)

    

# %%
idx = np.argsort(freqs0)
freqs0 =  freqs0[idx]
ps0 = ps0[idx]
ps_avg = []
f_avg = []
comp = 999999999
for ix, f in enumerate(freqs0):
    if f == comp:
        continue
    else:
        comp = f
        index = np.where(freqs0 == f)[0]
        ps_tmp = np.sum(ps0[index])/len(index)
        ps_avg.append(ps_tmp)
        f_avg.append(f)

plt.plot(f_avg, ps_avg)

# %%
