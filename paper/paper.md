---
title: 'DefectCorrectionsNotebook: Post-processing workflow for the application of finite-size corrections to periodic electronic structure calculations of charged defect supercells performed with the FHI-aims software package'
tags:
  - python3
  - point defects
  - electronic structure calculations
  - first principles
  - FHI-aims
  - charged defects
  - finite-size correction
authors:
  - name: Suzanne K. Wallace
    orcid: 0000-0003-4925-4768
    affiliation: "1,2"
  - name: Tong Zhu
    orcid: 
    affiliation: 3
  - name: Volker Blum
    orcid: 0000-0001-8660-7230
    affiliation: 3
affiliations:
 - name: Department of Chemistry, Centre for Sustainable Chemical Technologies, University of Bath, Claverton Down, Bath, BA2 7AY, UK
   index: 1
 - name: Department of Materials, Imperial College London, Exhibition Road, London SW7 2AZ, UK
   index: 2
 - name: Department of Mechanical Engineering and Materials Science, Duke University, Durham, North Carolina 27708, USA
   index: 3
date: 10 March 2019
bibliography: paper.bib
---

JOSS guidelines:
- A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience
- **Important: Note the paper should not include software documentation such as API (Application Programming Interface) functionality, as this should be outlined in the software documentation
- A clear statement of need that illustrates the purpose of the software
- Mentions (if applicable) of any ongoing research projects using the software or recent scholarly publications enabled by it


# Summary
The defect physics of a material often plays a decisive role in determining the performance of that material in various optoelectronic devices. First-principles electronic structure calculations of point defects can provide valuable insights into the defect physics of a material to complement experimental measurements [@Freysoldt2014]. The calculated formation energy, $\Delta H_{D,q}$ of a defect in a host material can be used to determine the expected equilibrium concentration of that type of defect when the material is synthesised under particular conditions. When a defect forms, there is a trade-off in the cost of breaking bonds and the gain in configurational entropy, $S$. The resulting crystal configuration will be that which minimises the Gibbs free energy, $G$, where $H$ is the enthalpy. 
$G = H - TS$
The trade-off between $H$ and $TS$ results in an equilibrium defect concentration, $n$, which depends on the formation energy of the defect and the number of lattice sites, $N$.
$n = N e^{\frac{-\Delta H_{D,q}}{k_B T}}$
$\Delta H_{D,q}$ directly depends on the chemical potential of the species involved in the creation of the defect, $\mu_i$. $\mu_i$ can be tuned experimentally as it is a function of temperature and pressure. Consequently, this allows for some tunability in the concentrations of particular defects in a material through the synthesis conditions. The formation energy of a defect can differ considerably when considering the defect in different charge states, it is therefore important to consider the different possible charge states when studying the defect physics of a material. Periodic boundary conditions (PBCs) are usually used by electronic software packages for computationally efficient simulations of solids. However, the use of PBCs can result in unphysical artifacts when simulating charged defects in the dilute (non-interacting) limit. As defects are often present in only parts per million of host lattice atoms, this is often the desired regime to simulate.


The supercell approach is a common method to overcome some of the limitations of PBCs for simulating dilute defects, by increasing the distance between defects and their periodic images. The formation energy, $\Delta H_{D,q=0}$, of charge neutral point defects in a supercell can be obtained by comparing the total energy calculated for the defective supercell to that of an equivalent perfect supercell of the host crystal and then considering the species added to or removed from the perfect supercell when the particular defect is formed. $\Delta H_{D,q=0}$ for a charge neutral defect is shown in the equation below where $E_{D,q=0}$ is the total energy of the defective supercell, $E_{host}$ is the total energy of an equivalent supercell of the perfect, bulk host crystal, $n_i$ is the number of atoms of species i added to ($n_i > 0 $) or removed from ($n_i < 0$) the chemical reservoir when the defect is formed. $\mu_i$ is the chemical potential of species i, referenced to the total energy of the pure element in its standard state, $E_i$. The chemical potential of a species i is the change in energy when one particle of type i is added to the system [@Zhang1991].
$ \Delta H_{D,q=0} = E_{D,q=0} - E_{host} + \sum_i n_i (E_i + \mu_i ) $
Additional complexities arise when attempting to obtain the formation energy of a charged isolated point defect. Firstly, there is a strong and long-ranged Coulomb interaction between charged supercells in PBCs and this converges slowly with increased supercell size.
Secondly, the charge of the defect system does not match that of the perfect bulk reference system. It is therefore necessary to introduce a chemical potential to account for the change in energy when electrons are added to or removed from the system when creating a defect in a given charge state.
Thirdly, electronic structure calculations with PBCs for a charged unit cell (effectively) include a neutralising homogeneous background charge to avoid infinite charge, which is not present in the calculation for the perfect equivalent supercell [@Freysoldt2014]. Consequently, the expression for the defect formation energy for a charge neutral defect, $\Delta H_{D,q=0}$, must be modified to that shown below for a charged defect in charge state q, $\Delta H_{D,q}$.
$ \Delta H_{D,q} = E_{D,q} - E_{host} + \sum_i n_i (E_i + \mu_i) + q[\epsilon_F + \epsilon_{\nu} + \Delta \nu_{0/b}] + E^q_{corr} $
Additional terms for a charged defect are: q (the charge state of the defect), $\epsilon_F$ (position of the Fermi level in the band gap),  $q \epsilon_{\nu}$ (energy of bulk VBM) and $\Delta \nu_{0/b}$ (term used to align the electrostatic potential of the valence band maximum for the bulk and defect supercells) and $E^q_{corr}$. The latter term usually represents multiple corrections performed after the periodic electronic structure calculations, one such correction is that for interactions between a charged defect and its periodic images, the 'image-charge' correction and another is the 'band filling' correction. 


- Give defect formation energy equation. Show neutral first then highlight additional terms/ corrections needed in charged case (see thesis). Cite several papers for more details on terms in this equation. (see thesis)
- Brief overview of first principles calculations of defects and finite-size errors for charged defects when using the supercell method with periodic DFT for use of well-established electronic structure software packages (see thesis - briefly mention iic, pa and bf corrections) + refer to literature for further reading (doi: 10.1103/RevModPhys.86.253)
- Highlight scalability of FHI-aims (cite) and importance of large defect supercells in many cases
- Mention that there are other existing workflows developed to perform finite-size correction schemes to outputs from other electronic structure softwares (cite pylada, pycdt, coffee) and that this workflow makes use of components of the CoFFEE python code (doi: 10.1016/j.cpc.2018.01.011), which was developed for outputs of the quantumESPRESSO software package (cite).


- This workflow has been designed to combine the benefits of step-by-step explanation and transparency of processing steps of a jupyter notebook with scriptability to allow for a convenient, automated application of finite-size corrections to a set of periodic electronic structure calculations for supercells containing charged point defects. This workflow has been developed to be compatible with outputs from the all-electron electronic structure software package FHI-aims (cite).
- The FNV post-processing scheme (doi: 10.1103/PhysRevLett.102.016402) is used to apply finite-size corrections to charged defect supercells when periodic electronic structure calculations have been performed with the FHI-aims software package.
- The potential alignment step (towards the end of the workflow) uses averaged atom-site potentials with the sampling region proposed by Kumagai and Oba (doi: 10.1103/PhysRevB.89.195205), see Fig. 2(a) in the publication. Currently, band filling corrections (most relevant for shallow defects with delocalised defect-induced charge) are not available in this workflow, only image-charge interaction and potential alignment corrections are applied. The addition of band filling corrections would be an example of a possible extension of this project.

# Defect processing example
The example below shows a visualisation generated with the VESTA (cite) software package for a supercell of ??? containing a ??? point defect (a) along with selected outputs from the notebook workflow, including plots of the defect coordinates located by the workflow (b) and alignment of the potentials of the defect and perfect host supercell (c).

**Show VESTA visual of defect supercell, then supercell location plot from notebook and final potential alignment plot? GaAs vacancy?**

# Acknowledgements

# References