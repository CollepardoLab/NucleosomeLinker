import numpy as np
from scipy.spatial.distance import cdist
import sys
from mymathfuncs import *
import copy

# linkers
LINKERS=np.loadtxt("linker_lengths.txt",dtype=int)

NF=len(LINKERS)-1


#Number of nucls to make
#NF=int(sys.argv[2])

# bp per bead
NB=5


DNA_FILE="AA_1kx5/ref_frames.dat"
CA_FILE="AA_1kx5/CAs.txt"

DCX=DCY="56"
DCZ="40"
DDX=DDY=DDZ="24"


# different ways to fit the NB base pairs per minimal bead
MY_QUAT_AV=2

def get_dna_bps():
    
    # get the dna data

    dna_data = list(open(DNA_FILE))

    # parse the data file
    # first line:
    # xxx base-pairs
    #
    # each bp organised like:
    # ...     1 G-C   # ...1>J:1169_:[DG5]G - ...1>K:1590_:[DC3]C\n'
    #   152.0681   242.6043   116.5436  # origin\n'
    #    -0.5248     0.6823    -0.5090  # x-axis\n'
    #     0.8288     0.5460    -0.1226  # y-axis\n'
    #     0.1943    -0.4862    -0.8520  # z-axis\n'

    n_dna = int(dna_data[0].split()[0])

    dna_beads = list()
    line = 1
    for n in range(n_dna):
        bps = dna_data[line].split()[2]
        line = line + 1
        coods = dna_data[line].split()
        r = np.array(
            [str(coods[0]),
             str(coods[1]),
             str(coods[2])],
            dtype=float)
        line = line + 1
        coods = dna_data[line].split()
        x = np.array(
            [str(coods[0]),
             str(coods[1]),
             str(coods[2])],
            dtype=float)
        line = line + 1
        coods = dna_data[line].split()
        y = np.array(
            [str(coods[0]),
             str(coods[1]),
             str(coods[2])],
            dtype=float)
        line = line + 1
        coods = dna_data[line].split()
        z = np.array(
            [str(coods[0]),
             str(coods[1]),
             str(coods[2])],
            dtype=float)
        line = line + 1

        # need to convert to a quaternion
        q = exyz_to_q(x, y, z)

        dna_beads.append([r, q])
    
    return dna_beads

def get_core():
    data=list(open(CA_FILE))
    
    coords=[]
    chains=[]
    for line in data:
        sline=line.split()
        x=np.array(sline[6:9],dtype=float)
        chain=sline[4]
        chains.append(chain)
        coords.append(x)

    return np.array(coords),chains

def fit_core(data):
    # first compute the reference frame vector then convert to quat

    # origin is the com
    com = get_com(data)

    # x axis points to specific beads
    x1=113
    x2=600
    x=(data[x1,:]+data[x2,:])*0.5
    xv = x-com
    
    #zef 
    z1=466
    z2=953

    zef1=data[z1,:]
    zef2=data[z2,:]

    zefv=zef1-zef2

    yv=-np.cross(xv,zefv)

    zv=np.cross(xv,yv)

    xv=unit_vec(xv)
    yv=unit_vec(yv)
    zv=unit_vec(zv)

    q=exyz_to_q(xv,yv,zv)

    return com,q


def trim(indna,L1,L2):
    dna=copy.deepcopy(indna)
    current = len(dna)
    print(current)
    
    nrl=147+L1+L2

    if nrl>current:
        #print("Error trying to have NRL of ",nrl ," which is greater than ",current," in the reference structure")
        #sys.exit()

        # need to add nrl:
        RISE=3.4
        TWIST=35*np.pi/180.0

        #left=-int(np.ceil((current-nrl)/2))
        #right=-int(np.floor((current-nrl)/2))
        left=L1
        right=L2

        print(left)
        print(right)

        # add to the start
        first_dna = dna[0]
        x = first_dna[0]
        q = first_dna[1]

        new_dna=list()

        for i in range(0,left):

            ex,ey,ez = q_to_exyz(q)

            new_x = x - ez*RISE
            new_q = quat_norm( quat_mul( quat_norm( quat_axis_angle(ez,-TWIST)),q))

            # need to get the positions of the phosphates
            #ex, ey, ez = q_to_exyz(new_q)
            #unit_vec(ey)
            #unit_vec(ez)

            #pos1 = (new_x - WIDTH * rotation(ey, ez, -P_ANGLE))

            #pos2 = (new_x + WIDTH * rotation(ey, ez, P_ANGLE))

            #new_dna.append(["P1", "_", pos1])
            #new_dna.append(["P2", "_", pos2])

            new_dna.append([new_x, new_q,"X-X"])

            q = quat_norm( new_dna[-1][1])
            x = new_dna[-1][0]


        new_dna.reverse()

        dna = new_dna + dna

        # add to the end

        last_dna = dna[-1]
        x = last_dna[0]
        q = quat_norm(last_dna[1])

        for i in range(0,right):

            ex,ey,ez = q_to_exyz(q)

            unit_vec(ez)

            new_x = x + ez*RISE
            new_q = quat_norm( quat_mul( quat_norm(quat_axis_angle(ez,TWIST)),q))

            #unit_vec(ey)
            dna.append([new_x, new_q,"X-X"])
            # need to get the positions of the phosphates
            #ex, ey, ez = q_to_exyz(quat_norm( dna_full[-1][-1]))

            #pos1 = (dna_full[-1][2] + WIDTH * rotation(ey, ez, P_ANGLE))

            #pos2 = (dna_full[-1][2] - WIDTH * rotation(ey, ez, -P_ANGLE))

            #dna_full.append(["P1", "_", pos1])
            #dna_full.append(["P2", "_", pos2])

            q = quat_norm(dna[-1][1])
            x = dna[-1][0]


        return dna



    else:    
    
        #left=int(np.ceil((current-nrl)/2))
        #right=int(np.floor((current-nrl)/2))
        left =L1
        right=L2

        print(current, left, right, current-left-right,nrl)

        if right>0:
            return dna[left:-right]
        else:
            return dna[left:]




#def trim(dna,nrl):
#    current = len(dna)
#
#    if nrl>current:
#        print("Error trying to have NRL of ",nrl ," which is greater than ",current," in the reference structure")
#        sys.exit()
#    
#    left=int(np.ceil((current-nrl)/2))
#    right=int(np.floor((current-nrl)/2))
#
#    print(current, left, right, current-left-right,nrl)
#
#    if right>0:
#        return dna[left:-right]
#    else:
#        return dna[left:]

def my_quat_av(qs):

    if qs.shape[0]==1:
        ret = np.mean(qs,axis=0)
        ret = ret/np.linalg.norm(ret)
        return ret

    elif MY_QUAT_AV==1:
        l=len(qs)
        i=l//2
        j=i-1
        
        q=qs[i]
        #q=slerp(qs[j],qs[i],[0.5])[0]
        ex,ey,ez=q_to_exyz(q)
        exx=ex
        # average the ez
        ezs=[]
        for q in qs:
            ex,ey,ez=q_to_exyz(q)
            ezs.append(ez)
        ezs=np.array(ezs)
        ez=np.mean(ezs,axis=0)
        ez=unit_vec(ez)
        ex=exx
        ey=np.cross(ez,ex)
        q=exyz_to_q(ex,ey,ez)
        return q
    else:

        qs=np.array(qs)
            
        q0=qs[0,:]
        for i in range(1,5):
            dd=np.dot(q0,qs[i,:])
            if dd<0: # double cover
                #print(q0,qs[i,:],"dd ",dd)
                qs[i,:]=-qs[i,:]
        
        ret = np.mean(qs,axis=0)
        ret = ret/np.linalg.norm(ret)
        #print(ret)
        return ret


def make_data(cores, all_dna,total_ids,total_types,total_molids):
    atoms=[]
    ells=[]
    new_bonds=[]
    bplines=[]
    total_xs=[]

    # first cores
    i=0
    for core in cores:
        com = core[0]
        total_xs.append(com)
        q = core[1]
        atoms.append(str(total_ids[i])+ " "+str(total_types[i])+" " + str(com[0])+" "+str(com[1])+" "+str(com[2])+" 1 1 "+str(total_molids[i])+" 0\n")
        ells.append(str(total_ids[i]) + " "+DCX+" "+DCY+" "+DCZ+" "+str(q[0])+" "+str(q[1])+" "+str(q[2])+" "+str(q[3])+"\n")
        i+=1

    # dna
    dids=[]
    for dna in all_dna:
        com=dna[0]
        total_xs.append(com)
        q=dna[1]
        atoms.append(str(total_ids[i])+ " "+str(total_types[i])+" " + str(com[0])+" "+str(com[1])+" "+str(com[2])+" 1 1 "+str(total_molids[i])+" 0\n")
        ells.append(str(total_ids[i]) + " "+DDX+" "+DDY+" "+DDZ+" "+str(q[0])+" "+str(q[1])+" "+str(q[2])+" "+str(q[3])+"\n")
        dids.append(total_ids[i])
        bplines.append(str(total_ids[i]) + " XX\n")
        i+=1
    
    # dna-dna bonds
    b=1
    for i in range(1,len(dids)):
        id1=dids[i-1]
        id2=dids[i]

        if not (total_types[id1-1] ==3 and total_types[id2-1]==3 and total_molids[id1-1]==total_molids[id2-1]):
        
            bond = [b,1,id1,id2]
            new_bonds.append(bond)
            b+=1

    bp_file = open("DNA_sequence.txt", "w")
    bp_file.write("# " + str(len(bplines)) + "\n")
    for line in bplines:
        bp_file.write(line)
    bp_file.close()


    total_xs=np.array(total_xs)
    xmin=np.min(total_xs[:,0])-500
    ymin=np.min(total_xs[:,1])-500
    zmin=np.min(total_xs[:,2])-500
    
    xmax=np.max(total_xs[:,0])+500
    ymax=np.max(total_xs[:,1])+500
    zmax=np.max(total_xs[:,2])+500
    
    bonds=[]
    for bond in new_bonds:
        bonds.append(" ".join(str(x) for x in bond)+"\n")

    nbonds=len(bonds)
    data_file = open("data.txt", "w")
    data_file.write("#LAMMPS data file\n")
    data_file.write("\n")
    data_file.write(str(len(atoms)) + " atoms\n")
    data_file.write(str(len(ells)) + " ellipsoids\n")
    data_file.write(str(nbonds) + " bonds\n")
    data_file.write("0 angles\n")
    data_file.write("0 dihedrals\n")
    data_file.write("0 impropers\n")
    data_file.write("\n")
    data_file.write("3 atom types\n")
    data_file.write("1 bond types\n")
    data_file.write("0 angle types\n")
    data_file.write("\n")
    data_file.write(str(xmin)+" "+str(xmax)+" xlo xhi\n")
    data_file.write(str(ymin)+" "+str(ymax)+" ylo yhi\n")
    data_file.write(str(zmin)+" "+str(zmax)+" zlo zhi\n")
    data_file.write("\n")
    data_file.write("Atoms\n")
    data_file.write("\n")
    for line in atoms:
        data_file.write(line)
    data_file.write("\n")
    data_file.write("Ellipsoids\n")
    data_file.write("\n")
    for line in ells:
        data_file.write(line)
    data_file.write("\n")
    data_file.write("Bonds\n")
    data_file.write("\n")
    for line in bonds:
        data_file.write(line)
    data_file.write("\n")
    
    data_file.close()

# def get_molids(dna,core,chains):
#     # strip tails from the core
#     tail_index= {"A":[[1, 40]],
#                 "B":[[1, 25]],
#                 "C":[[1, 20], [114, 128]],
#                 "D":[[1, 25]],
#                 "E":[[1, 40]],
#                 "F":[[1, 25]],
#                 "G":[[1, 20], [114, 128]],
#                 "H":[[1, 25]]}
#     tails=[]

#     k=1
#     for c in chains:
#         tail=False
#         idxs=tail_index[c]
#         #print(idxs)

#         for idx in idxs:
#             offset = chains.index(c)
#             #print(offset)
#             if k>=idx[0]+offset and k<=idx[1]+offset:
#                 tail=True 

#         tails.append(tail)
#         k+=1
    
#     tails = np.array(tails)
#     #print(tails)

#     core=core[~tails,:]
#     #print(core.shape)

#     # compute core-dna contacts
#     # need to include the phosphate locations

#     p1=[]
#     p2=[]
#     ds=[]
#     WIDTH = 8.5
#     P_ANGLE = 20.0*np.pi/180.0

#     for d in dna:
#         ex, ey, ez = q_to_exyz(d[1])

#         pos1 = (d[0] + WIDTH * rotation(ey, ez, P_ANGLE))

#         pos2 = (d[0] - WIDTH * rotation(ey, ez, -P_ANGLE))
#         p1.append(pos1)
#         p2.append(pos2)

#         ds.append(d[0])

#     p1=np.array(p1)
#     p2=np.array(p2)
#     ds=np.array(ds)

#     CUT=7.5

#     M1=cdist(p1,core)<CUT
#     M2=cdist(p2,core)<CUT
#     M3=cdist(ds,core)<CUT

#     M1=np.sum(M1,axis=1)
#     M2=np.sum(M2,axis=1)
#     M3=np.sum(M3,axis=1)

#     M=M1+M2+M3

#     idx=np.where(M>0)
#     print(idx)
#     i1 = idx[0][0]
#     i2 = idx[0][-1]
    
#     M[i1:i2+1]=1

#     print(M)

#     return M



if __name__ == "__main__":

    # get data from AA nucleosome 
    dna=get_dna_bps() # list of x,q pairs

    core,chains=get_core()


    # replace core atoms with the minimal bead
    ref_core_com,ref_core_q = fit_core(core)


    # program stages:
    # 1. trim to desired NRL still at bp resolution

    # 2. replicate nucleosomes to desired size

    # 3. remap to 5bp resolution using simple method


    # 1.

    ref_dna = dna


    trimmed_dna=trim(ref_dna,LINKERS[0],LINKERS[1])
    
    assert(len(trimmed_dna) == (147+LINKERS[0]+LINKERS[1]))


    #2. 

    
    all_dna=[*trimmed_dna]
    all_cores=[(ref_core_com,ref_core_q)]

    #left=-int(np.ceil((147-NRLS[0])/2))
    #right=-int(np.floor((147-NRLS[0])/2))
    left =LINKERS[0]
    right=LINKERS[1]
    
    print(left,right)
    dna_molid = [0 for i in range(left)] + [1 for i in range(147)] + [0 for i in range(right)]

    assert(len(dna_molid) == len(trimmed_dna))
    
    
    all_dna_molid=[*dna_molid]

    for i in range(1,NF):
        #nrl=NRLS[i]



        # transform all dna

        ## first need to trim it to the desired nrl

        trimmed_dna = trim(ref_dna,0,LINKERS[i+1])
        
        assert(len(trimmed_dna)==(147+0+LINKERS[i+1])) 
 


        # need to work out the transformation


        endx = all_dna[-1][0]
        endq = all_dna[-1][1]

        ex,ey,ez = q_to_exyz(endq)
        #twist = np.pi/180*float(open(NAFLEX).readlines()[3].split()[5]) 
        #rise = float(open(NAFLEX).readlines()[3].split()[2])
        twist=34.5*np.pi/180.0
        rise=3.5

        newez = ez
        newex = rotation(ex,ez,twist)
        newey = rotation(ey,ez,twist)

        newx = endx + rise*ez


        origr=trimmed_dna[0][0]
        origq=trimmed_dna[0][1]

        ex,ey,ez = q_to_exyz(origq)

        M1 = np.zeros((4,4))
        M1[:3,0]=ex
        M1[:3,1]=ey
        M1[:3,2]=ez
        M1[:3,3]=origr
        M1[3,3]=1

        M2 = np.zeros((4,4))
        M2[:3,0]=newex
        M2[:3,1]=newey
        M2[:3,2]=newez
        M2[:3,3]=newx
        M2[3,3]=1

        # M2 = T.M1
        # T = M2.inv(M1)

        T=np.matmul(M2,np.linalg.inv(M1))
        print(T)

        #left =-int(np.ceil((147-nrl)/2))
        #right=-int(np.floor((147-nrl)/2))
    
        left =0
        right=LINKERS[i+1]

        print(left,right)
        dna_molid = [0 for i in range(left)] + [1 for i in range(147)] + [0 for i in range(right)]

        assert(len(dna_molid) == len(trimmed_dna))
 
        for j in range(len(trimmed_dna)):
            M=np.zeros((4,4))
            r=trimmed_dna[j][0]
            q=trimmed_dna[j][1]
            ex,ey,ez=q_to_exyz(q)
            M[:3,0]=ex
            M[:3,1]=ey
            M[:3,2]=ez
            M[:3,3]=r
            M[3,3]=1

            new_M = np.matmul(T,M)
            new_r=new_M[:3,3]

            newex = new_M[:3,0]
            newey = new_M[:3,1]
            newez = new_M[:3,2]
            new_q = exyz_to_q(newex,newey,newez)

            all_dna.append((new_r,new_q))

            m=dna_molid[j]
            if m==1:
                m+=i

            all_dna_molid.append(m)

        # transform a core
        M=np.zeros((4,4))
        r=ref_core_com
        q=ref_core_q
        ex,ey,ez=q_to_exyz(q)
        M[:3,0]=ex
        M[:3,1]=ey
        M[:3,2]=ez
        M[:3,3]=r
        M[3,3]=1

        new_M = np.matmul(T,M)
        new_r=new_M[:3,3]

        newex = new_M[:3,0]
        newey = new_M[:3,1]
        newez = new_M[:3,2]
        new_q = exyz_to_q(newex,newey,newez)

        all_cores.append((new_r,new_q))


    

    # 3. remap resolution

    dna_types=[]
    new_dna=[]
    dna_molids=[]

    NM = len(all_dna)//NB

    print(NM)
    i=0
    for k in range(NM):
        print(k,NM)

        dna = all_dna[i:i+NB]

        xs = [d[0] for d in dna]
        qs = [d[1] for d in dna]

        xs=np.array(xs)
        qs=np.array(qs)

        x=get_com(xs)
        q=my_quat_av(qs)

        new_dna.append((x,q))

        # molids are 0 if linker or >0 if nucleosomal
        m=np.min(all_dna_molid[i:i+NB])
        if m>0:
            t=3
        else:
            t=2
        dna_types.append(t)
        dna_molids.append(m)

        i+=NB

    print(len(dna_types),len(new_dna),len(dna_molids))

    # correctly assign atom ttypes and ids
    ids = [i for i in range(1,len(new_dna)+len(all_cores)+1)]

    types = [1 for i in range(len(all_cores))] + dna_types
    molids = [i for i in range(1, len(all_cores)+1)] + dna_molids

    print(len(ids),len(types),len(molids))

    # 5. write the output file and make bond topology
    make_data(all_cores,new_dna,ids,types,molids)



