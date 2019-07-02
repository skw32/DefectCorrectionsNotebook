import numpy as np
import sys, string
import matplotlib.pyplot as plt    
import matplotlib   

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

def plot_atom_average_alignment(lattice,defect_charge,aims_atom_pots,model_pots,model,filename, user_ymin = None, user_ymax = None, user_xlabel = None, user_ylabel = None, user_title = None): 
    volume = np.dot(lattice[0,:],np.cross(lattice[1,:],lattice[2,:]))
    # Define the parameters corresponding to the Wigner-Seize Cells.   
    Rws = (volume*3.0/(4.0*np.pi))**(1.0/3.0)
    temp_a = (lattice[0][0]+lattice[1][0]+lattice[2][0])**2
    temp_b = (lattice[0][1]+lattice[1][1]+lattice[2][1])**2
    temp_c = (lattice[0][2]+lattice[1][2]+lattice[2][2])**2
    Rwl  = 0.5*np.sqrt(temp_a + temp_b + temp_c)
 
    plt.clf()
    matplotlib.rc('font',family='Times New Roman')
    plt.figure(figsize=(5.3,4.0))
  
    plt.xlabel(r'Distance from a defect ($\AA$)',fontsize=15,fontname = "Times New Roman")
    plt.ylabel('Potential (V)',fontsize=15,fontname = "Times New Roman")
    plt.xlim((0,Rwl))

    if model == False: 
        #Calculate the Potential Alignment Term
        Dim = len(aims_atom_pots[:,1])
        value = 0
        number = 0
        for i in range(Dim):
            dist = aims_atom_pots[i][0]
            if ((dist > Rws) and (dist < Rwl)):
                value += aims_atom_pots[i,1]
                number += 1.0

        value_sample = value/number
        value = -1.0*defect_charge*value/number
    
        value = '%.3f' % value
        pa_atom = value
 
        title = r'$\Delta E_{PA}^{LZ}(\alpha, q) = -q(V(\alpha,0)-V(Host))|_{far}$ = '+str(value) + ' eV'

        mean = (max(aims_atom_pots[:,1]) + min(aims_atom_pots[:,1]))/2.0

        # Options for user to override plot settings
        if user_title:
            title = user_title
        if user_xlabel:
            plt.xlabel(user_xlabel,fontsize=15,fontname = "Times New Roman")
        if user_ylabel:
            plt.ylabel(user_ylabel,fontsize=15,fontname = "Times New Roman")
        
        
        #Ploting the alignment plot       
        plt.axvline(x=Rws,linestyle='--',color='orange',lw =2.0) 
        plt.text(Rws,mean,r'$R_{ws}$',color='orange',fontsize=15,fontname = "Times New Roman")
    
        plt.title(title,fontsize=15,fontname = "Times New Roman")
    
        plt.plot(aims_atom_pots[:,0],aims_atom_pots[:,1],'o',ms=5,markerfacecolor='none', markeredgecolor='red',label=r'V($\alpha$,0)-V(Host)')
        plt.annotate(s='', xy=(Rwl,value_sample), xytext=(Rws,value_sample), arrowprops=dict(arrowstyle='<->',color = 'orange',lw = 2.0))
        plt.axhline(y=value_sample, xmin=Rws, xmax=Rwl, color='orange',label='Sampling Region')
        plt.legend()
        plt.savefig(filename,dpi=600)
        plt.show(block=False)
        return  pa_atom
    
    
    if model == True: 
        Diff = aims_atom_pots[:,1] - model_pots[:,1]
        Dim = len(Diff)
        number = 0
        value = 0

        for i in range(Dim):
            #print(Model[i,0])
            if ((model_pots[i,0] > Rws) and (model_pots[i,0] < Rwl)):
                value += Diff[i]
                number += 1.0

        value_sample = value/number 
        value = defect_charge*value/number 
        value = -1.0*value 
        value = '%.3f' % value
        pa_atom = value
        title = r'$-q \Delta $'+'='+str(value) 
        y_min = min(model_pots[:,1])
        y_max = max(model_pots[:,1])
        Ymin = y_min - 0.2*(y_max-y_min)
        Ymax = y_max + 0.2*(y_max-y_min)    
        mean = (Ymin+Ymax)/2.0

        # Options for user to override plot settings
        if user_ymin:
            Ymin = user_ymin
        if user_ymax:
            Ymax = user_ymax
        if user_title:
            title = user_title
        if user_xlabel:
            plt.xlabel(user_xlabel,fontsize=15,fontname = "Times New Roman")
        if user_ylabel:
            plt.ylabel(user_ylabel,fontsize=15,fontname = "Times New Roman")

        plt.ylim((Ymin,Ymax))    
        plt.axvline(x=Rws,linestyle='--',color='orange',lw =2.0) 
        plt.text(Rws,mean,r'$R_{ws}$',color='orange',fontsize=15,fontname = "Times New Roman")
        #plt.title(title,fontsize=15,fontname = "Times New Roman")
        plt.title('test 2',fontsize=15,fontname = "Times New Roman")
    
        #plt.plot(X1,Y1,'or',label='V(Defect,q)-V(Defect,0)')
        plt.plot(aims_atom_pots[:,0],aims_atom_pots[:,1],'o',ms=5,markerfacecolor='none', markeredgecolor='red',label=r'V($\alpha$,q)-V(Host)')
        plt.plot(model_pots[:,0],model_pots[:,1],'o',ms=5,markerfacecolor='none', markeredgecolor='blue',label=r'V(Model)')
        plt.plot(model_pots[:,0],Diff,'ok',ms=5,label=r'(V($\alpha$,q)-V(Host))-V(Model)')
        plt.annotate(s='', xy=(Rws,value_sample), xytext=(Rwl,value_sample), arrowprops=dict(arrowstyle='<->',color = 'orange',lw = 2.0))
        plt.axhline(y=value_sample, xmin=Rws, xmax=Rwl, color='orange',label='Sampling Region')
    
    
        plt.legend()
        plt.savefig(filename,dpi=600)
        
        return pa_atom  


def plot_planar_average_alignment(lattice_constant_z,Defect_pos,charge,X,aims_pots,model_pots,model,func,filename, user_ymin = None, user_ymax = None, user_xlabel = None, user_ylabel = None, user_title = None):
    # begin ploting planer_average_alignment  
    # find the postion which is most far away from the defect position.
    xmin = Defect_pos -0.5*lattice_constant_z
    xmax = Defect_pos + 0.5*lattice_constant_z
    if xmin <= 0 and xmax <= lattice_constant_z :
        corr_pos = xmax
    else:
        corr_pos = xmin
        
    idx = (np.abs(X-corr_pos)).argmin() 
    # this is the array index stands for the position which is most far away from the defect positions, we will use this to calculate the correction values lager .  

    #shift the defect position to zero in the potential alignment plot.
    dim = len(X)
    X_shift = np.zeros(dim)
    for i in range(dim):
        X_shift[i] = X[i] - Defect_pos
        if X_shift[i] <= 0:
            X_shift[i] += lattice_constant_z
    
    #define the appropriate ylimit  
    ymax = max(aims_pots)
    ymin = min(aims_pots) 
    if abs(ymax) < abs(ymin): 
        ymax = ymax
        ymin = -1.5*ymax
    else: 
        ymax = -1.5*ymin 
        ymin = ymin    
    y_min = ymin - 0.3*(ymax-ymin)
    y_max = ymax + 0.3*(ymax-ymin)

    if user_ymin:
        y_min = user_ymin
    if user_ymax:
        y_max = user_ymax

    #Generating the Plots    
    plt.clf()
    matplotlib.rc('font',family='Times New Roman')
    plt.figure(figsize=(5.3,4.0)) 
    plt.xlim((0,lattice_constant_z))
    plt.ylim((y_min,y_max)) 
    plt.xlabel(r'Z-coordinates ($\AA$)',fontsize=15,fontname = "Times New Roman")
    plt.ylabel('Potential (V)',fontsize=15,fontname = "Times New Roman")
    if user_xlabel:
        plt.xlabel(user_xlabel,fontsize=15,fontname = "Times New Roman")
    if user_ylabel:
        plt.ylabel(user_ylabel,fontsize=15,fontname = "Times New Roman")

    
    if model == False :  
        #Stands for the case which have not model potential contribution    
        value =  -1.0*charge*aims_pots[idx]
        value = '%.3f' % value 
        pa_planAv = value
        #print (value)
        title = r'$-q(V(\alpha,0)-V(Host))|_{far}$ = '+str(value) + ' eV'
 
        if user_title:
            title = user_title

        plt.title(title,fontsize=15,fontname = "Times New Roman")
    
        plt.axhline(y=aims_pots[idx], xmin=0, xmax=lattice_constant_z,linestyle = '--', lw=2,color='orange')
    
        plt.plot(X_shift,aims_pots,'or',ms=5,label=r'V($\alpha$,0)-V(Host)')
        plt.legend()
        plt.savefig(filename,dpi=600)
        plt.show(block=False)
        return pa_planAv
    else:  
        #Stands for the case which have the model potential contribution   
        b_A = 0.529177249    
        X_model = model_pots[:,0]*b_A
        Y_model = model_pots[:,1]

        #Move the model so that it actually at the defect positions. (fix the bug in the CoFFEE(numerical errors))   
        modelmax = max(abs(Y_model)) 
        dim_model = len(X_model)
        for i in range(dim_model): 
            if abs(Y_model[i]) == modelmax: 
                model_x_defect = X_model[i] 
        shift = model_x_defect - Defect_pos 
        X_model_temp = X_model - shift 
        for i in range(dim_model): 
            if X_model_temp[i] < 0: 
                X_model_temp[i] += lattice_constant_z 
            elif X_model_temp[i] > lattice_constant_z: 
                X_model_temp[i] -= lattice_constant_z 
    
        #sort the X and Y model so that it have the right sequence. 
        sorter = X_model_temp.argsort()
        Y_model_temp = Y_model[sorter]
        X_model_temp = np.sort(X_model_temp)  

#        print (Y_model_temp)   
        #generate a model potential at the exact same grid as that in FHI-aims output. (By polinomial fit, 12 dimension)  
        param = np.polyfit(X_model_temp, Y_model_temp, 12)
        poly = np.poly1d(param)
        dim = len(X)
        Y_model_corr = np.zeros(dim)
        for i in range(dim):
            Y_model_corr[i] = poly(X[i])
        
#        print (Y_model_corr)
        #Find the correction values 
        Diff = aims_pots - Y_model_corr
        value_sample = Diff[idx]
        value =  -1.0*charge*Diff[idx]
        value = '%.3f' % value 
        pa_planAv = value
        #   print (value)
        if func == 'CoFFEE': 
            title = r'-q((V($\alpha$,q)-V($\alpha$,0))-V(Model))|$_{far}$'+' = '+str(value) + ' eV'
            label1 = r'V($\alpha$,q)-V($\alpha$,0)'
            label2 = r'(V($\alpha$,q)-V($\alpha$,0))-V(Model)'
        if func == 'Frey': 
            title = r'-q((V($\alpha$,q)-V(Host))-V(Model))|$_{far}$'+' = '+str(value) + ' eV' 
            label1 = r'V($\alpha$,q)-V(Host)'
            label2 = r'(V($\alpha$,q)-V(Host))-V(Model)'
        # sorted all arrays   
        sorter_save = X_shift.argsort()
        aims_pots  = aims_pots[sorter_save]
        Y_model_corr = Y_model_corr[sorter_save]
        Diff = Diff[sorter_save]
        X_shift = np.sort(X_shift)
        
        if user_title:
            title = user_title

        plt.title(title,fontsize=15,fontname = "Times New Roman")
        
        plt.plot(X_shift,aims_pots,'or',label=label1)
        plt.plot(X_shift,Y_model_corr,'--b',label='V(Model)')
        plt.plot(X_shift,Diff,'--k',label = label2)
        plt.axhline(y=value_sample, xmin=0, xmax=lattice_constant_z,linestyle = '--',lw=2, color='orange')
        plt.legend()
        plt.savefig(filename,dpi=600)
        plt.show(block=False)
        return pa_planAv   




        # begin ploting planer_average_alignment without model potential contribution


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
