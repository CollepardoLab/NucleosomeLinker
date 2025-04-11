import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy.optimize as spo
import random
from scipy.interpolate import UnivariateSpline

NBINS=100
MINC=6.5735696879671400261e-02
MNUCL=200000
MDNA=650

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rcParams['font.family'] = 'sans-serif'
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)


# Function to convert CMYK to RGB
def cmyk_to_rgb(c, m, y, k):
    r = 255 * (1 - c) * (1 - k)
    g = 255 * (1 - m) * (1 - k)
    b = 255 * (1 - y) * (1 - k)
    return (r/255, g/255, b/255)

# Specific colors in CMYK and their lighter/darker versions
colors_cmyk = {
    "162": (0.40, 0.15, 0.00, 0.00),  # Lighter version of 172
    "167": (0.05, 0.45, 0.40, 0.00),  # Lighter version of 177
    "172": (0.80, 0.30, 0.00, 0.00),  # Original color
    "177": (0.10, 0.90, 0.80, 0.00),  # Original color
    "182": (0.90, 0.45, 0.10, 0.10)   # Darker version of 172
}

# Convert CMYK colors to RGB
colors_rgb = [cmyk_to_rgb(*cmyk) for nrl, cmyk in colors_cmyk.items()]


plt.figure(figsize=(5,4))
colors= ['#32b2ff', '#4cff32', '#9959ff', '#ff7f00','#f226f2','#e51932']
markers=["o","o","o","o","o","o"]
ls='--'

x_nrl=[]
y_nrl=[]

for nrl,C,m in zip(["172","173","174","175","176","177"],colors,markers+markers):
     
    VOLELEM=((120*120*800)/NBINS) / (MNUCL + (int(nrl)-147)*MDNA)
    
    # load data
    vals = np.loadtxt(nrl+"/data_points.txt")
    crit_point = np.loadtxt(nrl+"/crit_points.txt")
    curves = np.loadtxt(nrl+"/curves.txt")

    c1x = curves[:,0]/VOLELEM* 6.022e23 * (1e-24)
    c1y = curves[:,1]/MINC
    c2x = curves[:,2]/VOLELEM* 6.022e23 * (1e-24)
    c2y = curves[:,3]/MINC
    
    cd = crit_point[0]/VOLELEM* 6.022e23 * (1e-24)
    cds = crit_point[1]/VOLELEM* 6.022e23 * (1e-24)
    cc = crit_point[2]/MINC
    ccs = crit_point[3]/MINC

    plt.plot(c1x,c1y,linestyle=ls,color=C,linewidth=0.75)
    plt.plot(c2x,c2y,linestyle=ls,color=C,linewidth=0.75)
    plt.fill_betweenx(curves[:, 1] / MINC, c1x, c2x, color=C, alpha=0.3)    
    plt.errorbar(vals[:,1]/VOLELEM* 6.022e23 * (1e-24), vals[:,0]/MINC,xerr=vals[:,2]/VOLELEM* 6.022e23 * (1e-24),capsize=0.75,linestyle="",marker=m,color=C,label=nrl.rstrip("_ALT"),markersize=3,linewidth=0.5)
    plt.errorbar(vals[:,3]/VOLELEM* 6.022e23 * (1e-24), vals[:,0]/MINC,xerr=vals[:,4]/VOLELEM* 6.022e23 * (1e-24),capsize=0.75,linestyle="",marker=m,color=C,markersize=3,linewidth=0.5)
    plt.errorbar([cd],[cc],xerr=[cds],yerr=[ccs],capsize=0.75,linestyle="none",marker="x",mfc='none',color=C,markersize=3,linewidth=0.5)

    x_nrl.append(int(nrl))
    y_nrl.append(np.array((cc,ccs)))

plt.legend()
plt.xlabel(r"g/L")

plt.ylabel(r"$c/c_{crit}^{162}$")

plt.xlim((0,130))

plt.savefig('LLPS_fine.svg', format='svg')
plt.savefig("LLPS_nrl.pdf",bbox_inches="tight")

plt.tight_layout()
y_nrl=np.array(y_nrl)

plt.show()

