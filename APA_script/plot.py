# SUZY EDITS, PLAN:
## Take functions from the top as separate script
## Incorporate 'main' section below into notebook using notebook parameters as function args
## (Replace command line args for bottom of script with notebook parameters also)
## Generate plot in notebook and also save .png into defects_output_dir
#################################################

import matplotlib.pyplot as plt
import random
import time
from scipy import optimize
import numpy as np
import os
import sys
import math
import csv      

#main 
#charge = float(sys.argv[1])
#shift_H = float(sys.argv[2]) 
#shift_D = float(sys.argv[3])
#shift_H =  -15.02300189   
#shift_D =  -15.24318524 

### Suzy edits: put below into notebook to generate plot

### Suzy edits: Update these to take notebook params/ use dsa functions instead
defect_info, charge, shift_H, shift_D = ptt.sys_basic('system.in') # Suzy edits: replace this with notebook parameters??
atom_no,defect_no = geo_basic('geometry-host.in','geometry-defect.in')  
#print "atom_No, Defect_No, defect_line:", atom_no, defect_no 
defect_line,lattice_vec, host_corr = geo_compare(atom_no,defect_no,'geometry-host.in','geometry-defect.in')
result = apa.atomic_pot_fhiaims_plot(atom_no,defect_no,defect_line,lattice_vec,host_corr,'pot-host.out','pot-defect.out',shift_H,shift_D)
np.savetxt('atom_potential_FHI-aims',result)
volume = np.dot(lattice_vec[0,:],np.cross(lattice_vec[1,:],lattice_vec[2,:]))  


Rws = (volume*3.0/(4.0*np.pi))**(1.0/3.0)
temp_a = (lattice_vec[0][0]+lattice_vec[1][0]+lattice_vec[2][0])**2 
temp_b = (lattice_vec[0][1]+lattice_vec[1][1]+lattice_vec[2][1])**2 
temp_c = (lattice_vec[0][2]+lattice_vec[1][2]+lattice_vec[2][2])**2 
Rwl  = 0.5*np.sqrt(temp_a + temp_b + temp_c)
Model = np.loadtxt('atom_potential_model')
Diff = result[:,1] - Model[:,1] 
Dim = len(Diff)
number = 0
value = 0
for i in range(Dim): 
    if ((Model[i,0] > Rws) and (Model[i,0] < Rwl)): 
        value += Diff[i]
        number += 1.0 
value = charge*value/number 
value = -1.0*value 
value = '%.2f' % value 
title = defect_info[0] + ': '+r'$-\Delta q$'+'='+str(value) 
y_min = min(Model[:,1])
y_max = max(Model[:,1])
Ymin = y_min - 0.5*(y_max-y_min)
Ymax = y_max + 0.5*(y_max-y_min)      
plt.clf()
plt.xlim((0,Rwl))
plt.ylim((Ymin,Ymax))
plt.axvline(x=Rws,linestyle='--',color='orange') 
plt.title(title,fontsize=15,fontname = "Times New Roman")
#plt.plot(X1,Y1,'or',label='V(Defect,q)-V(Defect,0)')
plt.plot(result[:,0],result[:,1],'o',ms=5,markerfacecolor='none', markeredgecolor='red',label='V(Defect,q)-V(Host)')
plt.plot(Model[:,0],Model[:,1],'o',ms=5,markerfacecolor='none', markeredgecolor='blue',label='V(Model)')
plt.plot(Model[:,0],Diff,'ok',ms=5,label='(V(Defect,q)-V(Host))-V(Model)')
plt.legend()
plt.xlabel(r'Distance from Defect ($\AA$)',fontsize=15,fontname = "Times New Roman")
plt.ylabel('Potential (eV)',fontsize=15,fontname = "Times New Roman")
plt.savefig(FNV_atom_pa_file,dpi=600)