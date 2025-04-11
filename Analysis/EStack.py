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


# Twist deformation energy of DNA linker length
# from helical parameters in
# https://lcvmwww.epfl.ch/publications/data/articles/88/bj_85_2872.pdf


def Etwist_HEL(n_bp, kT):  # Defining the function

    # Diagonal stiffness constant in kcal/(mol*degree²)

    Ktw = np.mean(
        [0.0227, 0.0210, 0.0357, 0.0441, 0.0482, 0.0461, 0.0422, 0.0463, 0.0489, 0.0421]
    )
    # Estimate the twist deviation for a linker DNA of n_bp base pairs
    turns = n_bp / 10.2  # turns per linker DNA
    dev_turns_floor = np.mod(turns, 1)  # DeltaTurns wrt to lower integer
    dev_turns_top = 1 - np.mod(turns, 1)  # DeltaTurns wrt to upper integer
    # Minimum shift needed to recover flat nucleosomes
    dev_turns = np.min([dev_turns_floor, dev_turns_top])
    # Twist deviation in radians per basepair:
    # Delta_tw = dev_turns * pi / n_bp
    Delta_tw = dev_turns * 360 / n_bp
    # Toltal twist deformation energy in Kcal_mol per linker
    Etw_HP = 0.5 * Ktw * n_bp * (Delta_tw) ** 2
    return Etw_HP


def Etwist(n_bp, kT):  # Defining the function
    # DNA torsional persistence length in nm
    Lt = 75
    # linker DNA length in nm:
    l0 = 0.34
    # DNA torsional rigidity = Lt*kB*T/l0 in kcal/mol nm
    # s = 0.5 * Lt * kT / l0
    # Torsional rigidity of DNA in (kcal/mol):
    s = kT * Lt / (n_bp * l0)
    # Estimate the twist for each linker DNA of n_bp base pairs
    turns = n_bp / 10.4  # turns per linker DNA
    dev_turns_floor = np.mod(turns, 1)  # DeltaTurns wrt to lower integer
    dev_turns_top = 1 - np.mod(turns, 1)  # DeltaTurns wrt to upper integer
    # Minimum shift needed to recover flat nucleosomes
    dev_turns = np.min([dev_turns_floor, dev_turns_top])
    # Twist deviation in radians per basepair:
    Delta_tw_rad = dev_turns * 2 * pi
    # Delta_tw = dev_turns * 360
    # Delta_tw_rad = np.radians(Delta_tw)
    # Toltal twist deformation energy in Kcal_mol per linker
    Etw = 0.5 * s * (Delta_tw_rad) ** 2
    print(n_bp, Delta_tw_rad)
    return Etw


def Eval_and_plot_Energies():
    # Boltzmann constant k_B in kcal_mol/K
    kB = 1.987204259e-3
    # Temperature in K
    T = 300
    kT = kB * T

    # Linker DNA length in bp
    linkers = np.array([x for x in range(15, 70)])

    # Eff is the Face-to-face Inter-nucleosome interaction energy in kcal/mol
    # The value does not seem to change the probability much,
    # is that correct and emerging from the torsional rigidity dominating the
    # behaviour or did I make a mistake?
    Eff = 0

    # Estimate the torsional energy contribution:
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
