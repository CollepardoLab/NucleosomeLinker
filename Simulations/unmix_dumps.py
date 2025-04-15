import numpy as np
#import matplotlib.pyplot as plt
import sys


"""
Demixes the replica trajectories into the one at T0 and T63


"""


def get_traj_indexs(t,table):
    # return array mapping traj index into temperature order
    K=len(table)
    
    for k in range(1,K):
        #print(k)
        zu=table[k,0]
        zd=table[k-1,0]

        #print(t,zu,zd)

        if t >= zd and t < zu:
            row=table[k-1,1:]
            return [np.where(row == T)[0][0] for T in range(NUM_REP)]

    return None



# 1 load in the master log file

LOG_NAME    = "log.lammps"
N_REPLICAS  = 16
LOG_COLS    = 12
LOG_PATTERN = "log.temper_"
DUMP_PATTERN = "dna_temper_"
DOLOGS=False
DODUMPS=True

mlog=open(LOG_NAME,"r")
lines=list()
for line in mlog:
    lines.append(line)
mlog.close()

table=list()
for line in lines:
    sline=line.split()
    try:
        int(sline[0])
        if len(sline)==(N_REPLICAS+1):
            table.append(np.array(sline,dtype=int))
    except:
        pass

table=np.array(table)
#for line in table:
#    print(line)
print(table)

NUM_REP = N_REPLICAS

ABC=[str(i) for i in range(16)]

#ABC=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "AA", "AB", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AJ", "AK", "AL", "AM", "AN", "AO", "AP", "AQ", "AR", "AS", "AT", "AU", "AV", "AW", "AX", "AY", "AZ", "BA", "BB", "BC", "BD", "BE", "BF", "BG", "BH", "BI", "BJ", "BK", "BL"]


# 2 sort logs
if DOLOGS:

    print("reading in thermo logs")

    # read in the thermo logs

    thermos=list()
    for R in range(NUM_REP):
        print(R)
        thermo=list()
        lfile = open(LOG_PATTERN+ABC[R]+".txt","r")
        for line in lfile:
            try:
                a=np.array(line.split(),dtype=float)
                if len(a) == LOG_COLS:
                    thermo.append(a)
            except ValueError:
                pass
        thermo = np.array(thermo)
        #plt.plot(thermo[:,0], thermo[:,-1],label=ABC[R])
        thermos.append(thermo)
    #plt.legend()
    #plt.show()
    #plt.figure()
    print("sorting thermo logs")

    L=len(thermos[0])

    K=len(table)

    sorted_thermos=list()
    for R in range(NUM_REP):
        sorted_thermos.append(list())

    last_k=0
    for i in range(L):
        t=thermos[0][i,0]
        print(t)
        #if t < table[0,0]:
        #    for R in range(NUM_REP):
        #        sorted_thermos[R].append(thermos[R][i,:])
        #else:
        for k in range(last_k,K):
            #print(k)
            zu=table[k,0]
            zd=table[k-1,0]
            
            if t >= zd and t < zu:
                for T in range(NUM_REP):
                    X=np.where(table[k-1,1:] == T)[0][0]
                    #print(T,X,table[k-1,:])
                    # X is which replica is at temp 0
                    #print(T,X,i)
                    try:
                        sorted_thermos[T].append(thermos[X][i,:])
                    except IndexError:
                        pass
                last_k=k
                break

    print("thermos sorted")
    print("printing unmixed thermo logs")

    for st,n in zip(sorted_thermos,range(N_REPLICAS)):
        print(n)
        st=np.array(st)
        np.savetxt("log_equil_"+"T"+str(n)+".txt",st)


    
# 3 sort dump files
if DODUMPS:
    print("reading in dumps")
    # read in and sort in one go
    output_dumps=[]
    files=[]
    reading=[True]*NUM_REP
    for R in range(NUM_REP):
        dfilename=DUMP_PATTERN+ABC[R]+".dump"
        print(dfilename)
        dfile = open(dfilename,"r")
        files.append(dfile)


        out_file = open("equil_T" +str(R)+".dump","w")
        output_dumps.append(out_file)

    print(files)

    while any(reading):
        # loop over all files
        tsteps=[]
        frames=[]
        for dfile,n in zip(files,range(NUM_REP)):
            # get the next frame
            frame=[]
            # first line
            line1 = dfile.readline()
            if line1 != "ITEM: TIMESTEP\n":
                print("no frame on "+str(n))
                reading[n]==False
                continue
            else:
                reading[n]==True
            line2 = dfile.readline()
            tstep = int(line2)
            tsteps.append(tstep)
            #print(tstep)
            line3 = dfile.readline()
            if line3 != "ITEM: NUMBER OF ATOMS\n":
                break
            line4 = dfile.readline()
            N=int(line4)
            #print(N)
            line5=dfile.readline()
            line6=dfile.readline()
            line7=dfile.readline()
            line8=dfile.readline()
            line9=dfile.readline()

            frame.append(line1)
            frame.append(line2)
            frame.append(line3)
            frame.append(line4)
            frame.append(line5)
            frame.append(line6)
            frame.append(line7)
            frame.append(line8)
            frame.append(line9)

            for i in range(N):
                linei = dfile.readline()
                frame.append(linei)
            
            frames.append(frame)

        # have each frame at this timestep
        # reorder
        # look in the table
        if not (all(x==tsteps[0] for x in tsteps)):
            print("traj steps not equal")
            print(tsteps)
            break

        if len(tsteps)!= NUM_REP or len(frames)!=NUM_REP:
            print("not all frames read in")
            break

        traj_t = tsteps[0]
        print(traj_t)
        
        jmap = get_traj_indexs(traj_t,table)
        print(jmap)

        if jmap==None:
            continue

        print(jmap)

        # dump in T order
        #for T in range(NUM_REP):
        for T in [0]:
            outfile = output_dumps[T]

            frame = frames[jmap[T]]

            for line in frame:
                outfile.write(line)

        # next timestep

    for outfile in output_dumps:
        outfile.close()

    for dfile in files:
        dfile.close()


    print("finished")
