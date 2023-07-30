import sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.size"] = 24

def data_input(path_E_Rs, path_phi_0s):
    E_Rs = []
    with open(path_E_Rs) as f:
        for line in f.readlines():
            row = []
            toks = line.split(' ')
            for tok in toks:
                try:
                    num = float(tok)
                except ValueError as e:
                    #print(e, file=sys.stderr)
                    continue
                row.append(num)
            E_Rs.append(row)
        f.close()
    E_Rs = np.array(E_Rs)

    phi_0 = []
    with open(path_phi_0s) as f:
        for line in f.readlines():
            row = []
            toks = line.split(' ')
            for tok in toks:
                try:
                    num = float(tok)
                except ValueError as e:
                    #print(e, file=sys.stderr)
                    continue
                row.append(num)
            phi_0.append(row)
        f.close()
    phi_0 = np.array(phi_0)

    return E_Rs, phi_0

def make_plot_E_R(E_Rs):
    _, ax = plt.subplots(figsize = (12, 12))
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$E_{\mathrm{R}}$")
    ax.set_xlim(0, 50)
    ax.set_ylim(0, 1)
    ax.minorticks_on()

    ax.hlines([0.5], xmin=0, xmax=50,colors="gray")
    ax.scatter(E_Rs[:, 0], E_Rs[:, 1], c="black")

    plt.savefig('plot_E_R.jpg')

def make_plot_phi_0(phi_0):
    x = np.linspace(-5, 5, 101)
    y = np.pi**(-1 / 4) * np.exp(-x**2 / 2)

    _, ax = plt.subplots(figsize = (12, 12))
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$\phi_0(x)$")
    ax.set_xlim(-5, 5)
    ax.minorticks_on()

    ax.plot(x, y, c="gray")
    ax.scatter((phi_0[:, 0] + phi_0[:, 1]) / 2, phi_0[:, 2], c="black")

    plt.savefig('plot_phi_0.jpg')

path_E_Rs = "ER_evolution.txt"
path_phi_0s = "phi_0.txt"
E_Rs, phi_0 = data_input(path_E_Rs, path_phi_0s)

make_plot_E_R(E_Rs)
make_plot_phi_0(phi_0)