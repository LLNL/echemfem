import pandas
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr

df = pandas.read_csv("data.csv")

print(df)


def ideal_snes_slope(N): return df['num_processes'][0] * df['SNESSolve'][0] / N


def ideal_pcapply_slope(
    N): return df['num_processes'][0] * df['PCApply'][0] / N


fig = plt.figure()
ax = plt.axes(xscale='log', yscale='log')

ax.plot(df['num_processes'], df['SNESSolve'], marker='x')
ax.plot(df['num_processes'], [ideal_snes_slope(N)
                              for N in df['num_processes']], linestyle='dashed')

# plt.plot(df['num_processes'], df['PCApply'], marker='x')
# ax.plot(df['num_processes'], [ideal_pcapply_slope(N)
# for N in df['num_processes']], linestyle='dashed')

ax.set_yticks(df['SNESSolve'])
ax.get_yaxis().set_major_formatter(tkr.ScalarFormatter())

ax.set_xticks(df['num_processes'])
ax.get_xaxis().set_major_formatter(tkr.ScalarFormatter())

ax.minorticks_off()

ax.set_ylabel('solve time (s)')
ax.set_xlabel('number of processors')

fig.savefig('foo.png')
fig.clf()
