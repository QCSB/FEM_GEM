## GS GEM (modernized version) and microfluidics model

Dataset S1 (separate file). Modernized GS GEM provided in EXCEL format.  
Dataset S2 (separate file). Modernized GS GEM provided in SBML format.  
Dataset S3 (separate file). Modernized GS GEM provided in MAT format.  
Dataset S4 (separate file). Microfluidics full model for GS biofilm in Featool fea format.  
Dataset S5 (separate file). Microfluidics device structure in MAT format.  

GS GEM (modernized version): A modernized GS GEM has updated GEM for GS (GS_v2) (21) based on the original GS_v1 (22). First, any orphan metabolites and their associated reactions were removed from the model. Second, mass and charge balancing were performed by adding the metabolite information from the BiGG (Biochemical Genetic and Genomic) (23) database. If information about a metabolite could not be found, the relevant formula and charge were derived by mass and charge balancing of the reactions with which they are associated. For cytochromes, due to the absence of concrete knowledge, the formula and charge of the metabolites that were most extensively annotated in the BiGG database were used. Third, any erroneous names for the gene–protein-reaction (GPR) association were removed from its respective reaction, which was assigned an empty string for unknown GPR association instead of the reaction names. Finally, the updated GEM was converted to human-readable BiGG format (24). 

The modernized GS GEM structure is outlined in GS_v3.xlsx in Dataset S1. The new GEM of a GS cell (GS_v3) consists of 466 metabolites, 538 reactions, 889 genes and 367 reactions with gene–protein-reaction association. The improved GS GEM is provided in EXCEL (Dataset S1), SBML (Dataset S2), and MATLAB MAT (Dataset S3) format.

Model structure for the full microfluidics device (Featool and Matlab version): The full microfluidics model (GEM-FRTM) was built in FEATool Multiphysics (https://www.featool.com/) version 1.15 for Matlab. The Featool model for the full microfluidics device and the model structure readable by Matlab are provided in Dataset S4 and S5 (GS_microfluidics_model.fea and GS_microfluidics.mat).

21.	K. Zhuang et al., Genome-scale dynamic modeling of the competition between Rhodoferax and Geobacter in anoxic subsurface environments. ISME J 5, 305-316 (2011).
22.	R. Mahadevan et al., Characterization of metabolism in the Fe(III)-reducing organism Geobacter sulfurreducens by constraint-based modeling. Applied and Environmental Microbiology 72, 1558-1568 (2006).
23.	C. J. Norsigian et al., BiGG Models 2020: multi-strain genome-scale models and expansion across the phylogenetic tree. Nucleic Acids Res 48, D402-D406 (2020).
24.	Z. A. King et al., BiGG Models: A platform for integrating, standardizing and sharing genome-scale models. Nucleic Acids Res 44, D515-522 (2016).

