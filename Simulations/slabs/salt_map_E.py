import sys
import numpy as np
import scipy.interpolate

def salt_map(salt):
    # return the value of sigma we use,
    # and which TX to take the initial structures from

    salts=[0.15,
           0.115,
           0.082,
           0.07,
           0.06,
           0.052,
           0.05,
           0.042]

    Ts=[0,
        3,
        7,
        9,
        11,
        12,
        13,
        15,]
    
    sigma=[24,
        24.5,
        25,
        26,
        27,
        28,
        30,
        34]
    
    epsilon=[0.4,
            0.375,
            0.35,
            0.3,
            0.25,
            0.2,
            0.1,
            0.01,
            ]
    
    s_to_T = scipy.interpolate.interp1d(salts,Ts)
    s_to_sigma = scipy.interpolate.interp1d(salts,sigma)
    s_to_epsilon = scipy.interpolate.interp1d(salts,epsilon)

    T=s_to_T(salt)
    sigma=s_to_sigma(salt)
    E=s_to_epsilon(salt)
    
    T=int(np.round(T))
    #print(salt,T,sigma,E)

    return T,sigma,E


T,s,E = salt_map(float(sys.argv[1]))
#print("T=",T,"sigma=",s,"E=",E)
#print("variable E1 equal ",E)
#print("variable A equal ",s)
print(E)

