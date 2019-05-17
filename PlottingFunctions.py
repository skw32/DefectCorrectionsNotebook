import numpy as np
import sys, string


# Functions from Tong's plotting scripts -----------------------------------------------------------------------------------

def py_read(file):
    vol = np.load(file)
    grid = np.shape(vol)
    return grid,vol

def read_file(Filename):  
    # Reading the plane-average-file from FHI-aims output, 
    # and obtain the free-average-contribution, and the raw data 
    F = open(Filename, "r")
    res = F.readline().split()
    dim=len(res)
    average_free = float(res[dim-1])
    F.readline()
    F.readline()
    data = F.readlines()
    dim = len(data)
    result  = np.zeros([dim,2])
    i = 0 
    for line in data:
        words = line.split()
        result[i][0] = float(words[0])
        result[i][1] = float(words[1])
        i=i+1
    return average_free, result

# Functions from CoFFEE plotting script ('Examples/plot_fit.py' written by Mit Naik) --------------------------------------

def compute_fit(C_u,L2,L3,L4):
    '''
    Computes the fitting polynomial.
    C_u takes the 3 model energies
    L2, L3 and L4 are 1/\Omega^{1/3} for the corresponding cells.
    '''
    A = np.array([[ L2, L2**3, 1.], [  L3, L3**3, 1.] , [  L4, L4**3, 1. ] ])
    A_inv =  np.linalg.inv(A)
    X_u = np.dot(A_inv,C_u)
    return X_u

# Functions from CoFFEE plotting script ('PotentialAlignment/Utilities/plavg.py' written by Mit Naik) ---------------------

# For reading in_V file
def read_input(file_name):
    cell_dim = 1.0
    fp = open(file_name,'r')
    lines = fp.readlines()
    for il in range(len(lines)):
        if "file_name" in lines[il]:
            w = lines[il].split("=")
            if len(w) < 2 or len(w) > 3:
                print("ERROR while parsing input file: %s, line: %d"%(file_name,il))
                sys.exit()
            file_inp = w[1].split()[0]
        if "file_type" in lines[il]:
            w = lines[il].split("=")
            if len(w) < 2 or len(w) > 3:
                print("ERROR while parsing input file: %s, line: %d"%(file_name,il))
                sys.exit()
            file_type = w[1].split()[0]
        if "plt_dir" in lines[il]:
            w = lines[il].split("=")
            if len(w) < 2 or len(w) > 3:
                print("ERROR while parsing input file: %s, line: %d"%(file_name,il))
                sys.exit()
            plt_dir = w[1].split()[0]
        if "factor" in lines[il]:
            w = lines[il].split("=")
            if len(w) < 2 or len(w) > 3:
                print("ERROR while parsing input file: %s, line: %d"%(file_name,il))
                sys.exit()
            factor = w[1].split()[0]
        if "cell_dim" in lines[il]:
            w = lines[il].split("=")
            if len(w) < 2 or len(w) > 3:
                print("ERROR while parsing input file: %s, line: %d"%(file_name,il))
                sys.exit()
            cell_dim = eval(w[1])
    return file_inp,file_type,plt_dir,factor,cell_dim

def write2file(file_name,A,v_a):
    fp = open(file_name,'w')
    if len(A) != len(v_a):
        print("Error: len(A) != len(v_a)")
    for i in range(len(A)):
        fp.write("%4.3f %4.8f\n"%(A[i],v_a[i]))
    fp.close()

def pl_avg_a3(vol,a1_dim,a2_dim,a3_dim,step_l,factor, hartree, rydberg):
    A3 = []
    vol_a3 = np.zeros((a3_dim))
    for k in range(a3_dim):
        Sum1 = 0.
        for i in range(a1_dim):
            for j in range(a2_dim):
                Sum1 = Sum1 + vol[i][j][k]
        vol_a3[k] = Sum1/(a2_dim*a1_dim)
        A3.append(k*step_l)
    if factor == "Ryd":
        vol_a3 = vol_a3*rydberg
    elif factor == "Hartree":
        vol_a3 = vol_a3*hartree
    return vol_a3, np.array(A3)