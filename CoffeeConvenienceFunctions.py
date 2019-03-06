import subprocess
import os
from DefectSupercellAnalyses import *

def run_CoFFEE_solver(path_to_coffee_dir, defect_outputs_dir, super_x, super_y, super_z, defect_geom, sigma, cutoff, defect_charge, defect_x, defect_y, defect_z, dielectric_xx, dielectric_yy, dielectric_zz):
    write_CoFFEE_in_file(defect_outputs_dir, super_x, super_y, super_z, defect_geom, sigma, cutoff, defect_charge, defect_x, defect_y, defect_z, dielectric_xx, dielectric_yy, dielectric_zz)
    # Run CoFFEE Poisson solver to obtain E_q^per,m
    # Use subprocess.check_output to ensure coffee.py finishes running before notebook is allowed to proceed
    try:
        subprocess.check_output(str(os.path.join(path_to_coffee_dir, "coffee.py"))+" "+str(os.path.join(defect_outputs_dir,"cm_"+str(super_x)+"x"+str(super_y)+"x"+str(super_z)+"_in"))+" > "+str(os.path.join(defect_outputs_dir,"cm_"+str(super_x)+"x"+str(super_y)+"x"+str(super_z)+".out")), shell=True, stderr=subprocess.STDOUT)
    except Exception as err:
        logger.info(err.output)

def write_CoFFEE_in_V_file(super_z):
    in_V = open("in_V", "w")
    in_V.write("&plavg\n")
    in_V.write("file_name = V_r.npy\n")
    in_V.write("file_type = python\n")
    in_V.write("plt_dir = a3\n")
    in_V.write("factor =  None\n")
    in_V.write("cell_dim = "+str(super_z)+"\n") # z-dimension in units of Bohr
    in_V.write("/\n")
    in_V.close()

def write_CoFFEE_in_file(defect_outputs_dir, super_x, super_y, super_z, geom_file, sigma, cutoff, defect_charge, defect_x, defect_y, defect_z, dielectric_xx, dielectric_yy, dielectric_zz):
    ### Writing 'in' file for coffee.py
    coffee_in = open(os.path.join(defect_outputs_dir, "cm_"+str(super_x)+"x"+str(super_y)+"x"+str(super_z)+"_in"), "w")
    ## CELL PARAMS
    coffee_in.write("&CELL_PARAMETERS\n")
    coffee_in.write("\n")
    coffee_in.write("Lattice_Vectors(normalized):\n")
    # Read in lattice vectors from geometry.in, normalize and write to 'in' file for CoFFEE
    x_vecs, y_vecs, z_vecs = read_lattice_vectors(geom_file)
    a1_tot = x_vecs[0]+y_vecs[0]+z_vecs[0]
    a2_tot = x_vecs[1]+y_vecs[1]+z_vecs[1]
    a3_tot = x_vecs[2]+y_vecs[2]+z_vecs[2]
    coffee_in.write(str(x_vecs[0]/a1_tot)+"   "+str(y_vecs[0]/a1_tot)+"   "+str(z_vecs[0]/a1_tot)+"\n")
    coffee_in.write(str(x_vecs[1]/a2_tot)+"   "+str(y_vecs[1]/a2_tot)+"   "+str(z_vecs[1]/a2_tot)+"\n")
    coffee_in.write(str(x_vecs[2]/a3_tot)+"   "+str(y_vecs[2]/a3_tot)+"   "+str(z_vecs[2]/a3_tot)+"\n")
    coffee_in.write("\n")
    coffee_in.write("Cell_dimensions angstrom\n")
    coffee_in.write(str(a1_tot*super_x)+"  "+str(a2_tot*super_y)+"   "+str(a3_tot*super_z)+"\n")
    coffee_in.write("\n")
    coffee_in.write("Ecut="+str(cutoff)+" Hartree\n")
    coffee_in.write("/\n")
    coffee_in.write("\n")
    ## DIELECTRIC PARAMS
    coffee_in.write("&DIELECTRIC_PARAMETERS Bulk\n")
    coffee_in.write("Epsilon1_a1 = "+str(dielectric_xx)+"\n")
    coffee_in.write("Epsilon1_a2 = "+str(dielectric_yy)+"\n")
    coffee_in.write("Epsilon1_a3 = "+str(dielectric_zz)+"\n")
    coffee_in.write("/\n")
    coffee_in.write("\n")
    ## GAUSSIAN PARAMS (used for charge model)
    coffee_in.write("&GAUSSIAN_PARAMETERS:\n")
    coffee_in.write("Total_charge = "+str(defect_charge)+"\n")
    coffee_in.write("Sigma = "+str(sigma)+"\n")
    # Centre of Gaussian is set as defect location
    coffee_in.write("Centre_a1 = "+str(defect_x/a1_tot)+"\n") # CHECK THIS IS CORRECT! (in crystal units?!)
    coffee_in.write("Centre_a2 = "+str(defect_y/a2_tot)+"\n")
    coffee_in.write("Centre_a3 = "+str(defect_z/a3_tot)+"\n")
    coffee_in.write("/\n")
    coffee_in.close()