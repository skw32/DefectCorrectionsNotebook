import numpy as np 
cimport numpy as np 
cimport cython
DTYPE = np.complex128
ctypedef np.complex128_t DTYPE_t
@cython.boundscheck(False)
@cython.cdivision(True)



def atomic_3d(
             int defect_line,
             int atom_number,
             double beta,
             long[:] grid,
             double[:,:] latvec,
             double[:,:] host_corr,
             complex[:,:,:] V_G,
             double[:,:,:] G1, 
             double[:,:,:] G2,
             double[:,:,:] G3):
    cdef int i,j,k,l,jmax,kmax,lmax 
    cdef double[:] G = np.zeros((3))
    cdef double[:] R = np.zeros((3))
    cdef double[:] V_atomic = np.zeros((atom_number-1))
    cdef complex second 
    cdef double  third 
    print "TZ test calculations"
    jmax = (grid[0]-1)/2
    kmax = (grid[1]-1)/2
    lmax = (grid[2]-1)/2 
    for i in range(0,defect_line): 
        print i, " point"
        R[0] = host_corr[i,0]*latvec[0,0] + host_corr[i,1]*latvec[1,0] + host_corr[i,2]*latvec[2,0]
        R[1] = host_corr[i,0]*latvec[0,1] + host_corr[i,1]*latvec[1,1] + host_corr[i,2]*latvec[2,1]
        R[2] = host_corr[i,0]*latvec[0,2] + host_corr[i,1]*latvec[1,2] + host_corr[i,2]*latvec[2,2]
        for j in range(-jmax,jmax+1):
            for k in range(-kmax,kmax+1):
                for l in range(-lmax,lmax+1): 
                    G[0] = G1[j+jmax][k+kmax][l+lmax]
                    G[1] = G2[j+jmax][k+kmax][l+lmax]
                    G[2] = G3[j+jmax][k+kmax][l+lmax]
                    second = np.exp(1j*np.dot(G,R))
                    second = second*V_G[j+jmax][k+kmax][l+lmax]
                    third = np.exp(-0.5*np.dot(G,G)*beta*beta)
                    V_atomic[i] += np.real(second*third) 
    for i in range(defect_line+1,atom_number): 
        print i-1, " point"
        R[0] = host_corr[i,0]*latvec[0,0] + host_corr[i,1]*latvec[1,0] + host_corr[i,2]*latvec[2,0]
        R[1] = host_corr[i,0]*latvec[0,1] + host_corr[i,1]*latvec[1,1] + host_corr[i,2]*latvec[2,1]
        R[2] = host_corr[i,0]*latvec[0,2] + host_corr[i,1]*latvec[1,2] + host_corr[i,2]*latvec[2,2]
        for j in range(-jmax,jmax+1):
            for k in range(-kmax,kmax+1):
                for l in range(-lmax,lmax+1): 
                    G[0] = G1[j+jmax][k+kmax][l+lmax]
                    G[1] = G2[j+jmax][k+kmax][l+lmax]
                    G[2] = G3[j+jmax][k+kmax][l+lmax]
                    second = np.exp(1j*np.dot(G,R))
                    second = second*V_G[j+jmax][k+kmax][l+lmax]
                    third = np.exp(-0.5*np.dot(G,G)*beta*beta)
                    V_atomic[i-1] += np.real(second*third) 
 
    return np.array(V_atomic)
