---
title: 'DefectCorrectionsNotebook: Convenient and explanatory application of finite-size corrections to periodic electronic structure calculations of charged defect supercells with FHI-aims'
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
This workflow is designed to combine the benefits of step-by-step explanation and transparency of processing steps of a jupyter notebook with scriptability to allow for a convenient, automated application of finite-size corrections to a set of periodic electronic structure calculations for supercells containing charged point defects. 
- Highlight importance of defects for determining the performance of various devices (cite)
- Brief overview of first principles calculations of defects and finite-size errors for charged defects when using the supercell method with periodic DFT for use of well-established electronic structure software packages (see thesis - briefly mention iic, pa and bf corrections) + refer to literature for further reading (doi: 10.1103/RevModPhys.86.253)


Paragraph2:
- This python workflow enables the use of the FNV post-processing scheme (doi: 10.1103/PhysRevLett.102.016402) to apply finite-size corrections to charged defect supercells when performing periodic electronic structure calculations with the all-electron electronic structure software package FHI-aims (cite). Highlight scalability of FHI-aims (cite) and importance of large defect supercells in many cases.
- Mention that there are other existing workflows developed to perform such corrections to outputs from other electronic structure softwares (cite pylada, pycdt, coffee) and that this workflow makes use of components of the CoFFEE python code (doi: 10.1016/j.cpc.2018.01.011), which was developed for outputs of the quantumESPRESSO software package (cite).
- The potential alignment step (towards the end of the workflow) uses averaged atom-site potentials with the sampling region proposed by Kumagai and Oba (doi: 10.1103/PhysRevB.89.195205), see Fig. 2(a) in the publication. Currently, band filling corrections (most relevant for shallow defects with delocalised defect-induced charge) are not available in this workflow, only image-charge interaction and potential alignment corrections are applied. The addition of band filling corrections would be an example of a possible extension of this project.


# Example application
Maybe mention Tong's work for defect benchmarking between different implementations of electronic structure calculations as example of 'ongoing research projects using the software'??

Give defect formation energy equation. Show neutral first then highlight additional terms/ corrections needed in charged case (see thesis). Cite several papers for more details on terms in this equation.

Show VESTA visual of defect supercell, then supercell location plot from notebook and final potential alignment plot?

Give equation for defect concentrations based on defect formation energy - e.g. of important use of corrected formation energies of charged defects and highlight that formation energies can differ considerably when considering defects in different charge states, hence it is important to consider the different possible charge states when studying the defect physics of a material. 

# Acknowledgements

# References