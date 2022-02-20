%% Base GUI Microfluidics Model was created by FEATool Multiphysics Version 1.12.2, Build 20.05.122
%% Created with MATLAB 9.1.0.441655 (R2016b) PCWIN64.
%% Base GUI Model has been converted to m scripts and modified for M-model coupling and spatial-temporal current generation
%% microfluidcs model for phase 2 but can be extended to phase 1 and phase 3

Code: 

Geobacter.mat -- M-model for Geobacter
ion_mobility_full_precise_protons_v20.fea -- Basal GUI microfluidics model created by FEATool Multiphysics
M_model_microfluidics.m -- main code that was modified from the basal GUI microfluidics model for M-model coupling and data processing and plots
myfunGenome_acetate_warm_xy_phase2.m -- subroutine to obtain acetate oxidation rates from M-model
myfunGenome_growth_warm_xy_phase2.m -- subroutine to obtain growth rates from M-model
myfunGenome_h_warm_xy_phase2.m -- subroutine to obtain proton production rates from M-model
myfunGenome_hco3_warm_xy_phase2.m -- subroutine to obtain HCO3 production rates from M-model
myflow_compartment1.m -- subroutine to obtain flow rates in headspace compartment
myflow_compartment2.m -- subroutine to obtain flow rates in other compartments other than biofilm
