%% Base GUI Model created by FEATool Multiphysics Version 1.12.2, Build 20.05.122
%% Created with MATLAB 9.1.0.441655 (R2016b) PCWIN64.
%% Base GUI Model has been converted to m scripts and modification for M-model coupling and spatial-temporal current generation
%% This model includes phase 2 current generation
%% microfluidcs model for phase 2

clear all

global ace_points
global VOLUME
global model
global id_hco3
global id_h
global id_ace
global id_fe3
global growth
global Vmax
global Km
global Vhco3
global Vh
global fluxsize
global t_step
global flux_constraint

global flow_compartments2
global flow_compartments1

global GS_Currents_Fe3_phase2_time_series

global Acetate_phase2_time_series
global Acetate_phase2

global Fluxall

% mkdir SimDataPhase1
  mkdir MmodelSimDataPhase2

% mkdir SimDataPhase3

%% Loading Featool microfluidcs model developed from GUI ion_mobility_full_precise_protons_v20.fea
load('ion_mobility_full_precise_protons_20.fea','-mat')
%% set up M-model and warm start
changeCobraSolver('gurobi', 'LP');
load Geobacter.mat
model = addReaction(model,'EX_hco3(e)','reactionFormula','hco3[c] <=> hco3[e]')
model = changeRxnBounds(model,'EX_ac(e)',0,'u');
model = changeRxnBounds(model,'ATPM',0.45,'l'); 
model = changeRxnBounds(model,'ATPM',0.45,'u');
    
    id_ace = findRxnIDs(model, {'EX_ac(e)'});
    id_fe3 = findRxnIDs(model, {'FERCYT'}); % focytcc[c] + fe3[e]  -> ficytcc[c] + fe2[e]
    id_hco3 = findRxnIDs(model, {'EX_hco3(e)'});
    id_h = findRxnIDs(model, {'EX_h(e)'});
    
solution=solveCobraLP(model);
model.basis=solution.basis;
            
fluxsize=size(solution.full);
    
%% parameter setup 

     Vmax=3.7;
     Km=0.024e-06;
     
     flux_constraint=0.3;
    cphase2=0.14e-06;

    VOLUME=0.08*3*2; % mm^3
    dc_biofilm=1;
    dc_compartments=4;
    flow_compartments1=357;
    flow_compartments2=476;
    
    t_step=0.002;
    tf=0.202;
    
    Biomass=50e-06;


%% Grid setup
    ace_points = size(fea.grid.s(fea.grid.s==3),2);

%% Equation settings.
fea.phys.cd.dvar = { 'c' };
fea.phys.cd.sfun = { 'sflag1' };
fea.phys.cd.eqn.coef = { 'dts_cd', 'd_ts', 'Time scaling coefficient', { '1', '1', '1', '1', '1', '1' };
                         'd_cd', 'D', 'Diffusion coefficient', { dc_compartments, dc_compartments, '((y+0.08)/0.08)*0.4+((x+9.5)/9.5)*0.4', dc_compartments, dc_compartments, dc_compartments };
                         'u_cd', 'u', 'Convection velocity in x-direction', { 'cv_compartments1', 'cv_compartments1', '0', 'cv_compartments2', 'cv_compartments1', 'cv_compartments1' };
                         'v_cd', 'v', 'Convection velocity in y-direction', { '0', '0', '0', '0', '0', '0' };
                         'r_cd', 'R', 'Reaction rate', { '0', '0', 'Vace*c2', '0', '0', '0' };
                         'c0_cd', 'c_0', 'Initial condition for c', { cphase2, cphase2, cphase2, cphase2, cphase2, cphase2 } };
fea.phys.cd.eqn.seqn = 'dts_cd*c'' - d_cd*(cx_x + cy_y) + (u_cd*cx_t + v_cd*cy_t) - Z*c*rho_eps + Z*cx*Vx + Z*cy*Vy= r_cd';
fea.phys.cd.eqn.sdiff = ' - d_cd*(cx_x + cy_y)';
fea.phys.cd.eqn.sconv = '(u_cd*cx_t + v_cd*cy_t)';
fea.phys.cd.eqn.vars = { 'Concentration, c', 'c';
                         'Concentration gradient, c', { 'cx', 'cy' } };
fea.phys.cd.prop.isaxi = 0;
fea.phys.cd.prop.artstab.id = 0;
fea.phys.cd.prop.artstab.id_coef = 0.5;
fea.phys.cd.prop.artstab.sd = 0;
fea.phys.cd.prop.artstab.sd_coef = 0.25;
fea.phys.cd.prop.artstab.iupw = 1;
fea.phys.cd.prop.active = [ 1, 1, 1, 1, 1, 1 ];
fea.phys.cd.prop.intb = 0;
fea.phys.cd2.dvar = { 'c2' };
fea.phys.cd2.sfun = { 'sflag1' };
fea.phys.cd2.eqn.coef = { 'dts_cd2', 'd_ts', 'Time scaling coefficient', { '1', '1', '1', '1', '1', '1' };
                          'd_cd2', 'D', 'Diffusion coefficient', { '0', '0', '0', '0', '0', '0' };
                          'u_cd2', 'u', 'Convection velocity in x-direction', { '0', '0', '0', '0', '0', '0' };
                          'v_cd2', 'v', 'Convection velocity in y-direction', { '0', '0', '0', '0', '0', '0' };
                          'r_cd2', 'R', 'Reaction rate', { '0', '0', 'mu*c2', '0', '0', '0' };
                          'c20_cd2', 'c2_0', 'Initial condition for c2', { '0', '0', Biomass, '0', '0', '0' } };
fea.phys.cd2.eqn.seqn = 'dts_cd2*c2'' - d_cd2*(c2x_x + c2y_y) + (u_cd2*c2x_t + v_cd2*c2y_t) = r_cd2';
fea.phys.cd2.eqn.sdiff = ' - d_cd2*(c2x_x + c2y_y)';
fea.phys.cd2.eqn.sconv = '(u_cd2*c2x_t + v_cd2*c2y_t)';
fea.phys.cd2.eqn.vars = { 'Concentration, c2', 'c2';
                          'Concentration gradient, c2', { 'c2x', 'c2y' } };
fea.phys.cd2.prop.isaxi = 0;
fea.phys.cd2.prop.artstab.id = 0;
fea.phys.cd2.prop.artstab.id_coef = 0.5;
fea.phys.cd2.prop.artstab.sd = 0;
fea.phys.cd2.prop.artstab.sd_coef = 0.25;
fea.phys.cd2.prop.artstab.iupw = 1;
fea.phys.cd2.prop.active = [ 1, 1, 1, 1, 1, 1 ];
fea.phys.cd2.prop.intb = 0;
fea.phys.es.dvar = { 'V' };
fea.phys.es.sfun = { 'sflag1' };
fea.phys.es.eqn.coef = { 'eps_es', 'epsilon', 'Permittivity', { '1', '1', '1', '1', '1', '1' };
                         'Px_es', 'P_x', 'Polarization, x-component', { '0', '0', '0', '0', '0', '0' };
                         'Py_es', 'P_y', 'Polarization, y-component', { '0', '0', '0', '0', '0', '0' };
                         'rho_es', 'rho', 'Space charge density', { '0', '0', '0', '0', '0', '0' };
                         'V0_es', 'V_0', 'Initial condition for V', { '0', '0', '0', '0', '0', '0' } };
fea.phys.es.eqn.seqn = ' - eps_es*(Vx_x - Px_es_x + Vy_y - Py_es_y) = rho_es';
fea.phys.es.eqn.vars = { 'Electric potential, V', 'V';
                         'Electric field', 'sqrt(Vx^2+Vy^2)';
                         'Electric field, x-component', '-Vx';
                         'Electric field, y-component', '-Vy';
                         'Electric displacement', 'sqrt((-eps_es*Vx+Px_es)^2+(-eps_es*Vy+Py_es)^2)';
                         'Electric displacement, x-component', '-eps_es*Vx+Px_es';
                         'Electric displacement, y-component', '-eps_es*Vy+Py_es';
                         'Electric energy density', '0.5*((-eps_es*Vx+Px_es)*(-Vx)+(-eps_es*Vy+Py_es)*(-Vy))';
                         'Electric field, -grad(V)', { '-Vx', '-Vy' };
                         'Electric displacement', { '-eps_es*Vx+Px_es', '-eps_es*Vy+Py_es' } };
fea.phys.es.prop.isaxi = 0;
fea.phys.es.prop.active = [ 1, 1, 1, 1, 1, 1 ];
fea.phys.es.prop.intb = 0;
fea.phys.cd3.dvar = { 'c3' };
fea.phys.cd3.sfun = { 'sflag1' };
fea.phys.cd3.eqn.coef = { 'dts_cd3', 'd_ts', 'Time scaling coefficient', { '1', '1', '1', '1', '1', '1' };
                          'd_cd3', 'D', 'Diffusion coefficient', { dc_compartments, dc_compartments, '((y+0.08)/0.08)*0.4+((x+9.5)/9.5)*0.4', dc_compartments, dc_compartments, dc_compartments };
                          'u_cd3', 'u', 'Convection velocity in x-direction', { 'cv_compartments1', 'cv_compartments1', '0', 'cv_compartments2', 'cv_compartments1', 'cv_compartments1' };
                          'v_cd3', 'v', 'Convection velocity in y-direction', { '0', '0', '0', '0', '0', '0' };
                          'r_cd3', 'R', 'Reaction rate', { '0', '0', 'Vbicarbonate*c2', '0', '0', '0' };
                          'c30_cd3', 'c3_0', 'Initial condition for c3', { '0', '0', '0', '0', '0', '0' } };
fea.phys.cd3.eqn.seqn = 'dts_cd3*c3'' - d_cd3*(c3x_x + c3y_y) + (u_cd3*c3x_t + v_cd3*c3y_t) - Z*c*rho_eps + Z*c3x*Vx + Z*c3y*Vy = r_cd3';
fea.phys.cd3.eqn.sdiff = ' - d_cd3*(c3x_x + c3y_y)';
fea.phys.cd3.eqn.sconv = '(u_cd3*c3x_t + v_cd3*c3y_t)';
fea.phys.cd3.eqn.vars = { 'Concentration, c3', 'c3';
                          'Concentration gradient, c3', { 'c3x', 'c3y' } };
fea.phys.cd3.prop.isaxi = 0;
fea.phys.cd3.prop.artstab.id = 0;
fea.phys.cd3.prop.artstab.id_coef = 0.5;
fea.phys.cd3.prop.artstab.sd = 0;
fea.phys.cd3.prop.artstab.sd_coef = 0.25;
fea.phys.cd3.prop.artstab.iupw = 1;
fea.phys.cd3.prop.active = [ 1, 1, 1, 1, 1, 1 ];
fea.phys.cd3.prop.intb = 0;
fea.phys.cd4.dvar = { 'c4' };
fea.phys.cd4.sfun = { 'sflag1' };
fea.phys.cd4.eqn.coef = { 'dts_cd4', 'd_ts', 'Time scaling coefficient', { '1', '1', '1', '1', '1', '1' };
                          'd_cd4', 'D', 'Diffusion coefficient', { '4', '4', '((y+0.08)/0.08)*0.4+((x+9.5)/9.5)*0.4', '4', '4', '4' };
                          'u_cd4', 'u', 'Convection velocity in x-direction', { 'cv_compartments1', 'cv_compartments1', '0', 'cv_compartments2', 'cv_compartments1', 'cv_compartments1' };
                          'v_cd4', 'v', 'Convection velocity in y-direction', { '0', '0', '0', '0', '0', '0' };
                          'r_cd4', 'R', 'Reaction rate', { '0', '0', 'Vproton*c2', '0', '0', '0' };
                          'c40_cd4', 'c4_0', 'Initial condition for c4', { '0', '0', '0', '0', '0', '0' } };
fea.phys.cd4.eqn.seqn = 'dts_cd4*c4'' - d_cd4*(c4x_x + c4y_y) + (u_cd4*c4x_t + v_cd4*c4y_t) - Z*c*rho_eps + Z*c4x*Vx + Z*c4y*Vy  = r_cd4';
fea.phys.cd4.eqn.sdiff = ' - d_cd4*(c4x_x + c4y_y)';
fea.phys.cd4.eqn.sconv = '(u_cd4*c4x_t + v_cd4*c4y_t)';
fea.phys.cd4.eqn.vars = { 'Concentration, c4', 'c4';
                          'Concentration gradient, c4', { 'c4x', 'c4y' } };
fea.phys.cd4.prop.isaxi = 0;
fea.phys.cd4.prop.artstab.id = 0;
fea.phys.cd4.prop.artstab.id_coef = 0.5;
fea.phys.cd4.prop.artstab.sd = 0;
fea.phys.cd4.prop.artstab.sd_coef = 0.25;
fea.phys.cd4.prop.artstab.iupw = 1;
fea.phys.cd4.prop.active = [ 1, 1, 1, 1, 1, 1 ];
fea.phys.cd4.prop.intb = 0;
%% %% Constants and expressions.
fea.expr = {'cv_compartments1', 'myflow_compartment1(t)'; 'cv_compartments2', 'myflow_compartment2(t)';'Vace', 'myfunGenome_acetate_warm_xy_phase2(c,c2,x,y,t,cv_compartments1)';'mu', 'myfunGenome_growth_warm_xy_phase2(c,c2,x,y,t)'; 'Vbicarbonate', 'myfunGenome_hco3_warm_xy_phase2(c,c2,x,y,t)'; 'Vproton','myfunGenome_h_warm_xy_phase2(c,c2,x,y,t)';'Z', '0.17';'rho_eps', '0' };

%% Boundary settings.
fea.phys.cd.bdr.sel = [ 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 2, 3, -1, -1, -1, -1, -1, -1, -1 ];
fea.phys.cd.bdr.coef = { 'bcr_cd', 'c = c_0', 'Concentration', { 'c_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, cphase2, 0, 0, 0, 0, 0, 0, 0, cphase2, 0, 0, 0, 0, 0, 0, 0, 0 };
                         'bcc_cd', 'n.(-Dgrad c) = 0', 'Convective flux/outflow', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                         'bci_cd', 'n.(-Dgrad c + uc) = 0', 'Insulation/symmetry', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd+ny*v_cd)*c' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                         'bcf_cd', '-n.(-Dgrad c + uc) = N_0', 'Flux condition', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd+ny*v_cd)*c' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' } };
fea.phys.cd.bdr.coefi = { 'bcic_cd', 'n.(F_1-F_2) = 0, F_i=-D_igrad c_i  + u_ic_i', 'Continuity', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bcij_cd', '-n.(F_1-F_2) = N_0, F_i=-D_igrad c_i  + u_ic_i', 'Flux discontinuity', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd+ny*v_cd)*c' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                          'bcir_cd', 'c = c_0', 'Concentration', { 'c_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
fea.phys.cd.bdr.vars = { 'Normal diffusive flux of c', '-d_cd*(nx*cx+ny*cy)';
                         'Normal convective flux of c', '(nx*u_cd+ny*v_cd)*c';
                         'Normal total flux of c', '-d_cd*(nx*cx+ny*cy)+(nx*u_cd+ny*v_cd)*c' };
fea.phys.cd.prop.intb = 0;
fea.phys.cd2.bdr.sel = [ 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, -1, -1, -1, -1, -1, -1, -1 ];
fea.phys.cd2.bdr.coef = { 'bcr_cd2', 'c2 = c2_0', 'Concentration', { 'c2_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bcc_cd2', 'n.(-Dgrad c2) = 0', 'Convective flux/outflow', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bci_cd2', 'n.(-Dgrad c2 + uc2) = 0', 'Insulation/symmetry', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd2+ny*v_cd2)*c2' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                          'bcf_cd2', '-n.(-Dgrad c2 + uc2) = N_0', 'Flux condition', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd2+ny*v_cd2)*c2' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' } };
fea.phys.cd2.bdr.coefi = { 'bcic_cd2', 'n.(F_1-F_2) = 0, F_i=-D_igrad c2_i  + u_ic2_i', 'Continuity', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                           'bcij_cd2', '-n.(F_1-F_2) = N_0, F_i=-D_igrad c2_i  + u_ic2_i', 'Flux discontinuity', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd2+ny*v_cd2)*c2' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                           'bcir_cd2', 'c2 = c2_0', 'Concentration', { 'c2_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
fea.phys.cd2.bdr.vars = { 'Normal diffusive flux of c2', '-d_cd2*(nx*c2x+ny*c2y)';
                          'Normal convective flux of c2', '(nx*u_cd2+ny*v_cd2)*c2';
                          'Normal total flux of c2', '-d_cd2*(nx*c2x+ny*c2y)+(nx*u_cd2+ny*v_cd2)*c2' };
fea.phys.cd2.prop.intb = 0;
fea.phys.es.bdr.sel = [ 3, 3, 3, 3, 3, 1, 3, 3, 3, 1, 3, 3, -1, -1, -1, -1, -1, -1, -1 ];
fea.phys.es.bdr.coef = { 'bcp_es', 'V = V_0', 'Electric potential', { 'V_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { [], 0, '0', 0, 0, '250', 0, 0, 0, '0', 0, '0', 0, 0, 0, 0, 0, 0, 0 };
                         'bcg_es', 'V = 0', 'Ground/antisymmetry', [], { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                         'bcf_es', '-n.(-epsilongrad V+P) = rho_s', 'Surface charge', { 'rho_s' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '-nx*Px_es-ny*Py_es' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                         'bci_es', 'n.(-epsilongrad V+P) = 0', 'Insulation/symmetry', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)' } };
fea.phys.es.bdr.coefi = { 'bcic_es', 'n.(F_1-F_2) = 0, F_i=-epsilon_igrad V_i+P_i', 'Continuity', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bcij_es', '-n.(F_1-F_2) = rho_s, F_i=-epsilon_igrad V_i+P_i', 'Flux discontinuity', { 'rho_s' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bcir_es', 'V = V_0', 'Electric potential', { 'V_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
fea.phys.es.bdr.vars = { 'Electric potential', 'V';
                         'Surface charge', '-nx*(-eps_es*Vx+Px_es)-ny*(-eps_es*Vy+Py_es)' };
fea.phys.es.prop.intb = 0;
fea.phys.cd3.bdr.sel = [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, -1, -1, -1, -1, -1, -1, -1 ];
fea.phys.cd3.bdr.coef = { 'bcr_cd3', 'c3 = c3_0', 'Concentration', { 'c3_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bcc_cd3', 'n.(-Dgrad c3) = 0', 'Convective flux/outflow', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bci_cd3', 'n.(-Dgrad c3 + uc3) = 0', 'Insulation/symmetry', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd3+ny*v_cd3)*c3' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                          'bcf_cd3', '-n.(-Dgrad c3 + uc3) = N_0', 'Flux condition', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd3+ny*v_cd3)*c3' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' } };
fea.phys.cd3.bdr.coefi = { 'bcic_cd3', 'n.(F_1-F_2) = 0, F_i=-D_igrad c3_i  + u_ic3_i', 'Continuity', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                           'bcij_cd3', '-n.(F_1-F_2) = N_0, F_i=-D_igrad c3_i  + u_ic3_i', 'Flux discontinuity', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd3+ny*v_cd3)*c3' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                           'bcir_cd3', 'c3 = c3_0', 'Concentration', { 'c3_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
fea.phys.cd3.bdr.vars = { 'Normal diffusive flux of c3', '-d_cd3*(nx*c3x+ny*c3y)';
                          'Normal convective flux of c3', '(nx*u_cd3+ny*v_cd3)*c3';
                          'Normal total flux of c3', '-d_cd3*(nx*c3x+ny*c3y)+(nx*u_cd3+ny*v_cd3)*c3' };
fea.phys.cd3.prop.intb = 0;
fea.phys.cd4.bdr.sel = [ 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 2, 3, -1, -1, -1, -1, -1, -1, -1 ];
fea.phys.cd4.bdr.coef = { 'bcr_cd4', 'c4 = c4_0', 'Concentration', { 'c4_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, '0', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bcc_cd4', 'n.(-Dgrad c4) = 0', 'Convective flux/outflow', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bci_cd4', 'n.(-Dgrad c4 + uc4) = 0', 'Insulation/symmetry', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd4+ny*v_cd4)*c4' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                          'bcf_cd4', '-n.(-Dgrad c4 + uc4) = N_0', 'Flux condition', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd4+ny*v_cd4)*c4' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' } };
fea.phys.cd4.bdr.coefi = { 'bcic_cd4', 'n.(F_1-F_2) = 0, F_i=-D_igrad c4_i  + u_ic4_i', 'Continuity', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                           'bcij_cd4', '-n.(F_1-F_2) = N_0, F_i=-D_igrad c4_i  + u_ic4_i', 'Flux discontinuity', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd4+ny*v_cd4)*c4' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                           'bcir_cd4', 'c4 = c4_0', 'Concentration', { 'c4_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
fea.phys.cd4.bdr.vars = { 'Normal diffusive flux of c4', '-d_cd4*(nx*c4x+ny*c4y)';
                          'Normal convective flux of c4', '(nx*u_cd4+ny*v_cd4)*c4';
                          'Normal total flux of c4', '-d_cd4*(nx*c4x+ny*c4y)+(nx*u_cd4+ny*v_cd4)*c4' };
fea.phys.cd4.prop.intb = 0;
fea.phys.cd.bdr.sel = [ 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 2, 3, -1, -1, -1, -1, -1, -1, -1 ];
fea.phys.cd.bdr.coef = { 'bcr_cd', 'c = c_0', 'Concentration', { 'c_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, cphase2, 0, 0, 0, 0, 0, 0, 0, cphase2, 0, 0, 0, 0, 0, 0, 0, 0 };
                         'bcc_cd', 'n.(-Dgrad c) = 0', 'Convective flux/outflow', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                         'bci_cd', 'n.(-Dgrad c + uc) = 0', 'Insulation/symmetry', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd+ny*v_cd)*c' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                         'bcf_cd', '-n.(-Dgrad c + uc) = N_0', 'Flux condition', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd+ny*v_cd)*c' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' } };
fea.phys.cd.bdr.coefi = { 'bcic_cd', 'n.(F_1-F_2) = 0, F_i=-D_igrad c_i  + u_ic_i', 'Continuity', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bcij_cd', '-n.(F_1-F_2) = N_0, F_i=-D_igrad c_i  + u_ic_i', 'Flux discontinuity', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd+ny*v_cd)*c' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                          'bcir_cd', 'c = c_0', 'Concentration', { 'c_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
fea.phys.cd.bdr.vars = { 'Normal diffusive flux of c', '-d_cd*(nx*cx+ny*cy)';
                         'Normal convective flux of c', '(nx*u_cd+ny*v_cd)*c';
                         'Normal total flux of c', '-d_cd*(nx*cx+ny*cy)+(nx*u_cd+ny*v_cd)*c' };
fea.phys.cd.prop.intb = 0;
fea.phys.cd2.bdr.sel = [ 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, -1, -1, -1, -1, -1, -1, -1 ];
fea.phys.cd2.bdr.coef = { 'bcr_cd2', 'c2 = c2_0', 'Concentration', { 'c2_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bcc_cd2', 'n.(-Dgrad c2) = 0', 'Convective flux/outflow', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bci_cd2', 'n.(-Dgrad c2 + uc2) = 0', 'Insulation/symmetry', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd2+ny*v_cd2)*c2' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                          'bcf_cd2', '-n.(-Dgrad c2 + uc2) = N_0', 'Flux condition', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd2+ny*v_cd2)*c2' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' } };
fea.phys.cd2.bdr.coefi = { 'bcic_cd2', 'n.(F_1-F_2) = 0, F_i=-D_igrad c2_i  + u_ic2_i', 'Continuity', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                           'bcij_cd2', '-n.(F_1-F_2) = N_0, F_i=-D_igrad c2_i  + u_ic2_i', 'Flux discontinuity', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd2+ny*v_cd2)*c2' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                           'bcir_cd2', 'c2 = c2_0', 'Concentration', { 'c2_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
fea.phys.cd2.bdr.vars = { 'Normal diffusive flux of c2', '-d_cd2*(nx*c2x+ny*c2y)';
                          'Normal convective flux of c2', '(nx*u_cd2+ny*v_cd2)*c2';
                          'Normal total flux of c2', '-d_cd2*(nx*c2x+ny*c2y)+(nx*u_cd2+ny*v_cd2)*c2' };
fea.phys.cd2.prop.intb = 0;
fea.phys.es.bdr.sel = [ 3, 3, 3, 3, 3, 1, 3, 3, 3, 1, 3, 3, -1, -1, -1, -1, -1, -1, -1 ];
fea.phys.es.bdr.coef = { 'bcp_es', 'V = V_0', 'Electric potential', { 'V_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { [], 0, '0', 0, 0, '250', 0, 0, 0, '0', 0, '0', 0, 0, 0, 0, 0, 0, 0 };
                         'bcg_es', 'V = 0', 'Ground/antisymmetry', [], { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                         'bcf_es', '-n.(-epsilongrad V+P) = rho_s', 'Surface charge', { 'rho_s' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '-nx*Px_es-ny*Py_es' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                         'bci_es', 'n.(-epsilongrad V+P) = 0', 'Insulation/symmetry', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)', '(-nx*Px_es-ny*Py_es)' } };
fea.phys.es.bdr.coefi = { 'bcic_es', 'n.(F_1-F_2) = 0, F_i=-epsilon_igrad V_i+P_i', 'Continuity', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bcij_es', '-n.(F_1-F_2) = rho_s, F_i=-epsilon_igrad V_i+P_i', 'Flux discontinuity', { 'rho_s' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bcir_es', 'V = V_0', 'Electric potential', { 'V_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
fea.phys.es.bdr.vars = { 'Electric potential', 'V';
                         'Surface charge', '-nx*(-eps_es*Vx+Px_es)-ny*(-eps_es*Vy+Py_es)' };
fea.phys.es.prop.intb = 0;
fea.phys.cd3.bdr.sel = [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, -1, -1, -1, -1, -1, -1, -1 ];
fea.phys.cd3.bdr.coef = { 'bcr_cd3', 'c3 = c3_0', 'Concentration', { 'c3_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bcc_cd3', 'n.(-Dgrad c3) = 0', 'Convective flux/outflow', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bci_cd3', 'n.(-Dgrad c3 + uc3) = 0', 'Insulation/symmetry', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd3+ny*v_cd3)*c3' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                          'bcf_cd3', '-n.(-Dgrad c3 + uc3) = N_0', 'Flux condition', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd3+ny*v_cd3)*c3' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' } };
fea.phys.cd3.bdr.coefi = { 'bcic_cd3', 'n.(F_1-F_2) = 0, F_i=-D_igrad c3_i  + u_ic3_i', 'Continuity', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                           'bcij_cd3', '-n.(F_1-F_2) = N_0, F_i=-D_igrad c3_i  + u_ic3_i', 'Flux discontinuity', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd3+ny*v_cd3)*c3' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                           'bcir_cd3', 'c3 = c3_0', 'Concentration', { 'c3_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
fea.phys.cd3.bdr.vars = { 'Normal diffusive flux of c3', '-d_cd3*(nx*c3x+ny*c3y)';
                          'Normal convective flux of c3', '(nx*u_cd3+ny*v_cd3)*c3';
                          'Normal total flux of c3', '-d_cd3*(nx*c3x+ny*c3y)+(nx*u_cd3+ny*v_cd3)*c3' };
fea.phys.cd3.prop.intb = 0;
fea.phys.cd4.bdr.sel = [ 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 2, 3, -1, -1, -1, -1, -1, -1, -1 ];
fea.phys.cd4.bdr.coef = { 'bcr_cd4', 'c4 = c4_0', 'Concentration', { 'c4_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, '0', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bcc_cd4', 'n.(-Dgrad c4) = 0', 'Convective flux/outflow', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                          'bci_cd4', 'n.(-Dgrad c4 + uc4) = 0', 'Insulation/symmetry', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd4+ny*v_cd4)*c4' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                          'bcf_cd4', '-n.(-Dgrad c4 + uc4) = N_0', 'Flux condition', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd4+ny*v_cd4)*c4' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' } };
fea.phys.cd4.bdr.coefi = { 'bcic_cd4', 'n.(F_1-F_2) = 0, F_i=-D_igrad c4_i  + u_ic4_i', 'Continuity', [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
                           'bcij_cd4', '-n.(F_1-F_2) = N_0, F_i=-D_igrad c4_i  + u_ic4_i', 'Flux discontinuity', { 'N_0' }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { '(nx*u_cd4+ny*v_cd4)*c4' }, { '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0' };
                           'bcir_cd4', 'c4 = c4_0', 'Concentration', { 'c4_0' }, { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, [], { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
fea.phys.cd4.bdr.vars = { 'Normal diffusive flux of c4', '-d_cd4*(nx*c4x+ny*c4y)';
                          'Normal convective flux of c4', '(nx*u_cd4+ny*v_cd4)*c4';
                          'Normal total flux of c4', '-d_cd4*(nx*c4x+ny*c4y)+(nx*u_cd4+ny*v_cd4)*c4' };
fea.phys.cd4.prop.intb = 0;

%% Solver call.
fea = parsephys( fea );
fea = parseprob( fea );
disp('start ..........................................')
% fea.sol.settings.hook = @my_solve_hook % option to add electron neutrality constraint
%% Solver call.
tic
[fea.sol.u,fea.sol.t] = solvetime( fea, ...
                                   'iupw', [ 0, 0, 0, 0, 0 ], ...
                                   'linsolv', 'bicgstab', ...
                                   'tstep', t_step, ...
                                   'tmax', tf, ...
                                   'tstop', 1e-06, ...
                                   'ischeme', 2, ...
                                   'imass', 2, ...
                                   'icub', 'auto', ...
                                   'nlrlx', 1, ...
                                   'toldef', 1e-06, ...
                                   'tolchg', 1e-06, ...
                                   'reldef', 0, ...
                                   'relchg', 1, ...
                                   'maxnit', 20, ...
                                   'nproc', 4, ...
                                   'init', { 'c0_cd', 'c20_cd2', 'V0_es', 'c30_cd3', 'c40_cd4' }, ...
                                   'solcomp', [ 1, 1, 1, 1, 1, 1; 1, 1, 1, 1, 1, 1; 1, 1, 1, 1, 1, 1; 1, 1, 1, 1, 1, 1; 1, 1, 1, 1, 1, 1 ] );
toc
                    
%% current
% Plots for phase 1-3

time_phases2=0:t_step:tf;
Current_time_phases2=[GS_Currents_Fe3_phase2_time_series(1:end)];
result=[time_phases2' Current_time_phases2'];
csvwrite('plot_3.csv',result)
h=figure; plot(time_phases2,Current_time_phases2,'b-')
ylim([0 15])
xlim([0 0.20])
xlabel('Time (hr)'), ylabel('Current (uA)')


time_phases2=0:t_step:tf;
Ace_time_phases2=[Acetate_phase2_time_series(1:end)];
result=[time_phases2' Ace_time_phases2'];
h=figure; plot(time_phases2,Ace_time_phases2,'b-')
ylim([0 0.20e-06])
xlim([0 0.20])
xlabel('Time (hr)'), ylabel('Acetate (mmol mm^-^3)')

%% 

%% function for electron neutrality as an additional constraint

function [ A, f, prob ] = my_solve_hook( A, f, prob, solve_step )
%MY_SOLVE_HOOK Custom solver hook calling function.
%
%   [ A, F, PROB ] = SOLVE_HOOK( A, F, PROB, SOLVE_STEP )
%   Special solver hook function. Allows for special functions to be
%   called before and after the linear solver to modify the global
%   matrix A and right hand side/load vector F with the information in
%   the finite element problem struct PROB. SOLVE_STEP indicates if
%   the function call is before the linear solver call (=1), or to
%   clean-up after (=2).
%
%   A solve hook may be prescribed by entering a function name string
%   or handle in the fea.sol.settings.hook field. The corresponding
%   function will in each solve iteration then be called as
%
%       [A,f,fea] = fcn( A, f, fea, solve_step );

if( solve_step == 1 )
%  disp('Calling overall charge density to zero constraint function BEFORE interative solver call.');

  
  
  ind_c = find(strcmp(prob.dvar, 'c'));
    offset = sum(prob.eqn.ndof(1:ind_c-1));
  ind_dof_c = unique(prob.eqn.dofm{ind_c}(:)) + offset;
  n = length(ind_dof_c);
  
  ind_c3 = find(strcmp(prob.dvar, 'c3'));
  offset = sum(prob.eqn.ndof(1:ind_c3-1));
  ind_dof_c3 = unique(prob.eqn.dofm{ind_c3}(:)) + offset;
  
  ind_c4 = find(strcmp(prob.dvar, 'c4'));
  offset = sum(prob.eqn.ndof(1:ind_c4-1));
  ind_dof_c4 = unique(prob.eqn.dofm{ind_c4}(:)) + offset; 

  if( length(ind_dof_c) ~= length(ind_dof_c3) )
    error('c1 and c3 must have the same discretization (fea.sfun)')
  end

  if( length(ind_dof_c) ~= length(ind_dof_c4) )
    error('c1 and c4 must have the same discretization (fea.sfun)')
  end
  
  if( isstruct(A) )     % Sparsify matrix if necessary.
    A = sparse( A.rows, A.cols, A.vals, A.size(1), A.size(2) );
  end

  % Build constraint contributions.
  B = sparse(1, size(A,2));
  B(sub2ind(size(B), ones(1,n), ind_dof_c(:)')) = -1;
  B(sub2ind(size(B), ones(1,n), ind_dof_c3(:)')) = -1;
  B(sub2ind(size(B), ones(1,n), ind_dof_c4(:)')) = 1;
  % Add row to add constraint sum_i(c1_i + c2_i) = 0.
  A = [A, B';
       B, 0 ];
  f = [f; 0];
end

if( solve_step == 2 )
%  disp('Calling solver hook AFTER linear solver call.')

  % Remove system modifications.
  n0 = sum(prob.eqn.ndof);
  A = A(1:n0,1:n0);
  f = f(1:n0);

end

end