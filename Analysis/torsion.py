from pprint import PrettyPrinter

# Importing libraries
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from dataclasses import dataclass

# Set font for labels in coming figures
font16 = {
    "family": "serif",
    "color": "black",
    "weight": "normal",
    "size": 16,
}

# Set font for labels in coming figures
font20 = {
    "family": "serif",
    "color": "black",
    "weight": "normal",
    "size": 20,
}

# Twist deformation energy of DNA linker length from helical parameters


def Etwist_HEL(n_bp, kT):
    Ktw = np.mean(
        [0.0227, 0.0210, 0.0357, 0.0441, 0.0482, 0.0461, 0.0422, 0.0463, 0.0489, 0.0421]
    )
    turns = n_bp / 10.2
    dev_turns_floor = np.mod(turns, 1)
    dev_turns_top = 1 - np.mod(turns, 1)
    dev_turns = np.min([dev_turns_floor, dev_turns_top])
    Delta_tw = dev_turns * 360 / n_bp
    Etw_HP = 0.5 * Ktw * n_bp * (Delta_tw) ** 2
    
    # Print energy cost for each point
    print(f"Linker length: {n_bp} bp, Energy cost: {Etw_HP} kcal/mol (HEL method)")
    
    return Etw_HP


def Etwist(n_bp, kT):
    Lt = 75
    l0 = 0.34
    s = kT * Lt / (n_bp * l0)
    turns = n_bp / 10.4
    dev_turns_floor = np.mod(turns, 1)
    dev_turns_top = 1 - np.mod(turns, 1)
    dev_turns = np.min([dev_turns_floor, dev_turns_top])
    Delta_tw_rad = dev_turns * 2 * pi
    Etw = 0.5 * s * (Delta_tw_rad) ** 2
    
    # Print energy cost for each point
    print(f"Linker length: {n_bp} bp, Energy cost: {Etw} kcal/mol")
    
    return Etw


def Eval_and_plot_Energies():
    # Boltzmann constant k_B in kcal_mol/K
    kB = 1.987204259e-3
    T = 300
    kT = kB * T

    linkers = np.array([x for x in range(15, 70)])

    Eff = -4.8

    Et = np.array([Etwist(e, kT) for e in linkers])

    Et_HP = np.array([Etwist_HEL(e, kT) for e in linkers])

    Etotal_tw = np.add(Eff, Et)
    Etotal_twhp = np.add(Eff, Et_HP)
    PEtotal_tw = np.exp(-Etotal_tw / kT)
    PEtotal_twhp = np.exp(-Etotal_twhp / kT)

    Z_tw = np.max(PEtotal_tw)
    Z_twhp = np.max(PEtotal_twhp)

    plt.clf()

    plt.title("Face-to-face stacking", fontdict=font20)
    plt.xlabel("Linker Length (bp)", fontdict=font16)
    plt.ylabel("Normalised probability", fontdict=font16)

    plt.plot(
        linkers,
        PEtotal_tw / Z_tw,
        label="Exp. torsional persistence length",
        color="r",
        marker="o",
        linestyle="solid",
    )
    plt.plot(
        linkers,
        PEtotal_twhp / Z_twhp,
        label="Helical Parameters",
        color="k",
        marker="o",
        linestyle="solid",
    )
    plt.legend(fontsize=16)
    plt.show()

    plt.clf()

    plt.title(
        "Energy cost of twisting DNA to force face-to-face stacking", fontdict=font20
    )
    plt.xlabel("Linker Length (bp)", fontdict=font16)
    plt.ylabel("Energy (kcal/mol)", fontdict=font16)
    plt.plot(
        linkers,
        Etotal_tw,
        label="Exp. torsional persistence length",
        color="r",
        marker="o",
        linestyle="solid",
    )
    plt.plot(
        linkers,
        Etotal_twhp,
        label="Helical Parameters",
        color="k",
        marker="o",
        linestyle="solid",
    )
    plt.legend(fontsize=16)
    plt.savefig('Twisting.svg')
    plt.show()

    plt.clf()


Eval_and_plot_Energies()

