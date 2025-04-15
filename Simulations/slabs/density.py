import numpy as np
import sys
#import more_itertools
#import matplotlib.pyplot as plt

#Number of frames before equilibriation
#EQUIL=500#200
#CHUNKS=1
NBINS=100
ZMAX=8000
#NUM_AFRAMES = 1500


EQUILLIM=100

first=0.25
last=1.0
#LASTN=10

def get_z(fname):
    frames = list()
    traj_file = open(fname)
    I = 1
    frame = list()

    for line in traj_file:
        if line == "ITEM: TIMESTEP\n":
            if I%10==0:
                print(I)

            if I > 1:
                frame = np.array(frame)
                frames.append(frame)
                frame = list()
            I = I+1
            # if I>10:
            #     break

        else:  # try and read in line
            try:
                array_line = np.array(line.split(), dtype=float)
                if len(array_line) >3:
                    frame.append(array_line[3])
            except:
                pass

    traj_file.close()
    # last frame
    frame = np.array(frame)
    frames.append(frame)
    

    return frames


def wrap(z, zsize):
    a=zsize[0]
    b=zsize[1]
    d=b-a
    for i in range(len(z)):
        if z[i] > b:
            while z[i] > b:
                z[i]-=d
        if z[i]<a:
            while z[i] < a:
                z[i]+=d
    return z

frames=get_z(sys.argv[1])

zsize=(0,ZMAX)
print(len(frames))
C=(zsize[1]-zsize[0])/2


L=len(frames)
#TODO: How to decide which frames to use for analysis
#NUM_AFRAMES = 3*L//4
#LIM=EQUIL
#print(L,L-LIM)

#Ignore the last frame
#TODO: Temp changing to just looking at last few frames
#frames=frames[-1*(NUM_AFRAMES+1):-1]

#frames=frames[(len(frames)-LASTN):]

print(len(frames))

if len(frames)<EQUILLIM:
        print("not enough frames")
        sys.exit()

frames=frames[int(len(frames)*first):int(len(frames)*last)]

k=0

#chunks = [ list(x) for x in more_itertools.divide(CHUNKS,frames) ]
#n=0
#for frames in chunks:
all_zs=[]
#print(n,len(frames))
for zs in frames:
    # peridoic boundaries
    zs=wrap(zs,zsize)
    #need to center
    #print(k,len(frames)-L)
    k+=1
    # find mean, periodic, use angles
    # map boxl to 0-2pi
    zangs = zs*2*np.pi/zsize[1]
    meanangz = np.arctan2(np.mean(np.sin(zangs)),np.mean(np.cos(zangs)))
    #print(meanangz)
    meanz = meanangz*zsize[1]/(2*np.pi)
    #print(meanz)
    delta=C-meanz
    zs=zs+delta
    zs=wrap(zs,zsize)
    y,b=np.histogram(zs,bins=NBINS,range=zsize)
    x=(b[1:]+b[:-1])/2
    all_zs.append(y)
    #plt.plot(x,y)


all_zs=np.array(all_zs)
z=np.mean(all_zs,axis=0)
z_std = np.std(all_zs,axis=0)
#plt.figure()
#plt.plot(x,z,'rx-')

np.savetxt("dprof.txt",np.array([x,z,z_std]).T)

#n+=1
#plt.savefig("density.png")
    #plt.show()
