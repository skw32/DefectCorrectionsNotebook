"""
The following functions have been adapted with permission from pylada-defects
See: https://github.com/pylada/pylada-defects
and doi: 10.1016/j.commatsci.2016.12.040
"""

from pylada.crystal import Structure
from quantities import elementary_charge, eV, pi, angstrom
from pylada.physics import Ry, a0
from pylada.ewald import ewald
from pylada.crystal.defects import third_order

##############################################################
def get_madelungenergy(latt_vec_array: tuple, charge: int, epsilon: float, cutoff: float) -> float:
    """ Function returns leading first order correction term, i.e.,
        screened Madelung-like lattice energy of point charge
    Reference: M. Leslie and M. J. Gillan, J. Phys. C: Solid State Phys. 18 (1985) 973

    Args
        defect: pylada.vasp.Extract object
        charge: charge of point defect. Default 1e0 elementary charge
        epsilon: dimensionless relative permittivity, SKW: isotropic average of dielectric constant
        cutoff: Ewald cutoff parameter

    Returns
        Madelung (electrostatic) energy in eV                                                                                                                

    Note:
        1. Units in this function are either handled by the module Quantities, or\
        defaults to Angstrom and elementary charges
        2. Function is adopted from Haowei Peng's version in pylada.defects modules
    """

    ewald_cutoff = cutoff * Ry
    
    cell_scale = 1.0 # SKW: In notebook workflow cell parameters are converted to Cartesians and units of Angstroms
    # SKW: Create point charge in pylada.crystal.structure class (used for charge model)
    # http://pylada.github.io/pylada/userguide/crystal.html
    struc = Structure()
    struc.cell = latt_vec_array
    struc.scale = cell_scale
    struc.add_atom(0., 0., 0., "P", charge=charge)
    
    #Anuj_05/22/18: added "cutoff" in ewald syntax
    result = ewald(struc, cutoff=ewald_cutoff).energy / epsilon
    return -1*result.rescale(eV)

##############################################################
def thirdO(latt_vec_array: tuple, charge: int, n: int) -> float:
    """ Function returns 3rd order image charge correction, same as LZ fortran script
    Reference: S. Lany and A. Zunger, Phys. Rev. B 78, 235104 (2008)
               S. Lany and A. Zunger, Model. Simul. Mater. Sci. Eng. 17, 0842002 (2009)
               [Eq. 6, 7]

    Args
        defect: pylada.vasp.Extract object
        charge: charge of point defect. Default 1e0 elementary charge
        n: precision in integral of Eq. 7 (LZ 2009), larger the better

    Returns
        third image correction in eV
    """

    cell_scale = 1.0 # SKW: In notebook workflow cell parameters are converted to Cartesians and units of Angstroms  
    cell = (latt_vec_array*cell_scale) * angstrom.rescale(a0) 

    #Anuj_05/22/18:modified to "third_order"
    thirdO = third_order(cell, n) * (4e0*pi/3e0) * Ry.rescale(eV) * charge * charge

    return thirdO

##############################################################
def get_imagecharge(latt_vec_array: tuple, charge: int, epsilon: float, cutoff: float, n: int =20, verbose=True, **kwargs) -> tuple:
    """ Function returns complete image charge correction (Madelung + scaled 3rd Order)
        Reference: S. Lany and A. Zunger, Model. Simul. Mater. Sci. Eng. 17, 0842002 (2009)[Eq. 11]

    Args
        defect: pylada.vasp.Extract object
        charge: charge of point defect. Default 1e0 elementary charge
        epsilon: dimensionless relative permittivity, SKW: isotropic average of dielectric constant
        cutoff: Ewald cutoff parameter
        n: precision in integral of Eq. 7 (LZ 2009), larger the better
        verbose: True or False

    Returns
        Non-verbose = Madelung + scaled 3rd order image charge correction in eV
        Verbose = Madelung_energy, 3rd Order, shape-factor csh, scaling f, final_image_correction in eV
    """

    E1 = get_madelungenergy(latt_vec_array, charge=1e0, epsilon=1e0, cutoff=cutoff)
    E3 = -1.*thirdO(latt_vec_array, charge=1e0, n=n)

    if epsilon == 1e0:
        # epsilon==1e0, meaning vacuum                                                                                 
        print("epsilon=1e0, defect in vacuum. really!!")
        return E1/epsilon
    else:
        scaled_E3 = E3*(1e0 - 1e0/epsilon)
        csh = E3/E1
        f = csh*(1e0 - 1e0/epsilon)
        
        E_ic = (E1 + scaled_E3) * charge * charge /epsilon
        
        if verbose == False:
            return E_ic
        else:
            return ["{:0.3f}".format(float(E1)), "{:0.3f}".format(float(E3)), "{:0.3f}".format(float(csh)) \
                        ,"{:0.3f}".format(float(f)), "{:0.3f}".format(float(E_ic))]
