import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy.optimize as spo



ldllim=[500, 7500]
hdllim=[3800,4200]


dirs = glob.glob("./0.*")

dirs.sort()
print(dirs)

def fitfunc1(c,d,Cc):
    return d*(1.0-c/Cc)

def fitfunc2(x,b,Pc):
    return Pc+b*x


def get_crit_point(vals,errors):
    
    try:



        topn=int(np.loadtxt("topN.txt"))
        
        print(topn)

        print(vals)
    
        vals=vals[:topn,:]
        errors=errors[:topn,:]
    
    except OSError:
        pass

    # dont use points where vapor density is zero
    #valid  = vals[:,1]>0.0
    #
    
    #vals   = vals[valid,:]
    #errors = errors[valid,:]
    #vals=np.sort(vals,axis=0)
    
    #vals=vals[:4,:]

    print(vals)

    

    vf = vals[:,1]
    vf_error = errors[:,1]
    bf = vals[:,3]
    bf_error = errors[:,3]
    cf = vals[:,0]

    delta = bf-vf

    delta_error = np.sqrt((bf_error**2 + vf_error**2))
    deltap3 = delta**3.06
    deltap3_error = delta_error**3.06


    #popt,pcov = spo.curve_fit(fitfunc1,cf,deltap3, sigma=deltap3_error, absolute_sigma=True, bounds=((-1e15,0.06),(0,0.08)))
    popt,pcov = spo.curve_fit(fitfunc1,cf,deltap3, bounds=((-1e15,0.06),(0,0.08)))

    plt.figure()
    #plt.scatter(cf,deltap3)
    plt.errorbar(cf,deltap3,yerr=deltap3_error,capsize=5)
    plt.plot(cf,fitfunc1(cf, *popt))
    print(popt)

    crit_c = popt[1]
    d = popt[0]
    crit_c_std = np.sqrt(pcov[1,1])

    # crit d

    # critical d
    dsum = (bf+vf)/2.0

    dsum_error = np.sqrt(bf_error**2+vf_error**2)/2.0
    y=dsum



    x=crit_c-cf



    print(popt)
    plt.figure()
    plt.scatter(x,y)

    #popt,pcov = spo.curve_fit(fitfunc2,x,y,sigma=dsum_error,absolute_sigma=True,bounds=((-1e6,0),(0,1000)))
    popt,pcov = spo.curve_fit(fitfunc2,x,y,bounds=((-1e6,0),(0,1000)))

    print(popt)
    plt.plot(x,fitfunc2(x,*popt))

    crit_d = popt[1]
    s=popt[0]
    crit_d_std = np.sqrt(pcov[1,1])


    print(crit_c,crit_c_std,crit_d,crit_d_std)
        

    plt.figure()
    cx = np.linspace(crit_c,np.max(cf),1000)

    A = fitfunc1(cx,d,crit_c)**(1.0/3.06)
    B = fitfunc2((crit_c-cx),s,crit_d)

    rhoh = (A+2.0*B)/2.0
    
    rhol = rhoh - A

    plt.plot(rhoh,cx)
    plt.plot(rhol,cx)
    

    plt.scatter(bf,cf)
    plt.scatter(vf,cf)

    plt.scatter([crit_d],[crit_c])
    #plt.show()

    #plt.show()
    lines=[[rhoh,cx],[rhol,cx]]

    if len(bf) == 2:
            crit_c_std = 0.001
    return crit_c,crit_c_std,crit_d,crit_d_std,lines


vals=[]
errors=[]
for salt_dir in dirs:
        repeats=[]
        for rep in [1,2,3]:
                try:
                    dprof = np.loadtxt(salt_dir+"/"+str(rep)+"/dprof.txt")
                
                except OSError:
                    continue

                xx=dprof[:,0]
                d=dprof[:,1]

                plt.plot(xx,d,label=salt_dir)

                # get ldl and hld
                valid_ldl = (xx < ldllim[0])|(xx > ldllim[1])

                ldl_d = np.mean(d[valid_ldl])
                ldl_std = np.std(d[valid_ldl])
                
                valid_hdl = (xx > hdllim[0])&(xx < hdllim[1])

                hdl_d = np.mean(d[valid_hdl])
                hdl_std = np.std(d[valid_hdl])

                salt = float(salt_dir.lstrip("./"))

                val = (salt, ldl_d, ldl_std, hdl_d, hdl_std)
                print(val)
                repeats.append(val)
        
        repeats=np.array(repeats)
        
        print(repeats.shape)

        if repeats.ndim < 2:
                continue

        means = np.mean(repeats,axis=0)
        stds  = np.std(repeats,axis=0)

        print(means.shape)
        vals.append(means)
        errors.append(stds)

plt.axvline(x=ldllim[0])
plt.axvline(x=ldllim[1])

plt.axvline(x=hdllim[0])
plt.axvline(x=hdllim[1])




plt.legend()
#plt.show()





vals=np.array(vals)
errors=np.array(errors)

cc, ccs, cd, cds,lines = get_crit_point(vals,errors)


plt.figure()
#plt.plot(vals[:,1],vals[:,0])
#plt.plot(vals[:,3],vals[:,0])

plt.errorbar(vals[:,1], vals[:,0],xerr=errors[:,1],capsize=5,linestyle="none")
plt.errorbar(vals[:,3], vals[:,0],xerr=errors[:,3],capsize=5,linestyle="none")
plt.errorbar([cd],[cc],xerr=[cds],yerr=[ccs],capsize=5,linestyle="none",marker='*')
plt.plot(lines[0][0],lines[0][1],'k--')
plt.plot(lines[1][0],lines[1][1],'k--')

plt.gca().invert_yaxis()


# save the data
# data points
np.savetxt("data_points.txt",vals,header="conc, vapor, vapor_err,liquid, liquid_error")


# crit points
np.savetxt("crit_points.txt",np.array([cd,cds,cc,ccs]),header="crit_d, crit_d_err, crit_c, crit_c_err")

curve1=lines[0]
curve2=lines[1]
# curves
np.savetxt("curves.txt",np.array([curve1[0],curve1[1],curve2[0],curve2[1]]).T, header="c1x, c1y, c2x, c2y")


plt.savefig("PD.png")
#plt.show()


