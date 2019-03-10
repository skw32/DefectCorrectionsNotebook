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
The defect physics of a material often plays a decisive role in determining the performance of that material in various optoelectronic devices. First-principles electronic structure calculations of point defects can provide valuable insights into the defect physics of materials to complement experimental measurements [@Freysoldt2014]. The calculated formation energy, $\Delta H_{D,q}$ of a defect in a host material can be used to determine the expected equilibrium concentration of that type of defect when the material is synthesised under particular conditions. When a defect forms, there is a trade-off in the cost of breaking bonds and the gain in configurational entropy, $S$. The resulting crystal configuration will be that which minimises the Gibbs free energy, $G$, shown in Eq.~\ref{freeE} where $H$ is the enthalpy. 
$G = H - TS$
The trade-off between $H$ and $TS$ results in an equilibrium defect concentration, $n$, which depends on the formation energy of the defect and the number of lattice sites, $N$, as shown in Eq.~\ref{defect_conc2}.
\begin{equation}\label{defect_conc2}
n = N e^{\frac{-\Delta H_{D,q}}{k_B T}}
\end{equation}
$\Delta H_{D,q}$ directly depends on the chemical potential of the species involved in the creation of the defect, $\mu_i$. $\mu_i$ can be tuned experimentally as it is a function of temperature and pressure. Consequently, this allows for some tunability in the concentrations of particular defects in a material through the synthesis conditions. The formation energy of a defect can differ considerably when considering the defect in different charge states, it is therefore important to consider the different possible charge states when studying the defect physics of a material. Periodic boundary conditions are usually used by electronic software packages for computationally efficient simulations of solids. However, the use of periodic boundary conditions can result in unphysical artifacts when simulating charged defects in the dilute (non-interacting) limit.

The formation energy of a neutral defect is given by:
$$ $$
Where .. is .. and .. is .. The formation energy for a charged defect, however, is given by:
$$ $$
Additional terms in the expression for a charged defect are ...

- Give defect formation energy equation. Show neutral first then highlight additional terms/ corrections needed in charged case (see thesis). Cite several papers for more details on terms in this equation. (see thesis)
- Brief overview of first principles calculations of defects and finite-size errors for charged defects when using the supercell method with periodic DFT for use of well-established electronic structure software packages (see thesis - briefly mention iic, pa and bf corrections) + refer to literature for further reading (doi: 10.1103/RevModPhys.86.253)
- Highlight scalability of FHI-aims (cite) and importance of large defect supercells in many cases
- Mention that there are other existing workflows developed to perform finite-size correction schemes to outputs from other electronic structure softwares (cite pylada, pycdt, coffee) and that this workflow makes use of components of the CoFFEE python code (doi: 10.1016/j.cpc.2018.01.011), which was developed for outputs of the quantumESPRESSO software package (cite).

- This workflow has been designed to combine the benefits of step-by-step explanation and transparency of processing steps of a jupyter notebook with scriptability to allow for a convenient, automated application of finite-size corrections to a set of periodic electronic structure calculations for supercells containing charged point defects. This workflow has been developed to be compatible with outputs from the all-electron electronic structure software package FHI-aims (cite).
- The FNV post-processing scheme (doi: 10.1103/PhysRevLett.102.016402) is used to apply finite-size corrections to charged defect supercells when periodic electronic structure calculations have been performed with the FHI-aims software package.
- The potential alignment step (towards the end of the workflow) uses averaged atom-site potentials with the sampling region proposed by Kumagai and Oba (doi: 10.1103/PhysRevB.89.195205), see Fig. 2(a) in the publication. Currently, band filling corrections (most relevant for shallow defects with delocalised defect-induced charge) are not available in this workflow, only image-charge interaction and potential alignment corrections are applied. The addition of band filling corrections would be an example of a possible extension of this project.

# Example application
The example below shows a visualisation generated with the VESTA (cite) software package for a supercell of ??? containing a ??? point defect (a) along with selected outputs from the notebook workflow, including plots of the defect coordinates located by the workflow (b) and alignment of the potentials of the defect and perfect host supercell (c).

**Show VESTA visual of defect supercell, then supercell location plot from notebook and final potential alignment plot? GaAs vacancy?**

# Acknowledgements

# References