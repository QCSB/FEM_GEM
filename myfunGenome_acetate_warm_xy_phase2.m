function [ Vace ] = myfunGenome_acetate_warm_xy_phase2(c,c2,x,y,t,cv_compartments1)
% t;
global model
global id_ace
global id_fe3
global id_hco3
global id_h
global ace_points
global Vmax
global Km
global growth
global fluxsize
global deltat
global deltarace
global VOLUME
global Vhco3
global Vh
global t_step
global flux_constraint
global Fluxall
global GS_Currents_Fe3_phase2_time_series
global Acetate_phase2_time_series
global Acetate_phase2


c(c<0) = 0; 

Vmax=3.7;
Km=0.024e-06;

%% to test vmax change with time (optional)

 if (t(1) > 0.05) && (t(1) < 0.10)
    Vmax = 0.2+x/2+y*2;
  end

if t(1) > t_step
n=int16(t(1)/t_step+1);
n=n(1);
elseif t(1)==t_step
n=2;
else
n=1;
end


Vmax_ac_Numerator = (Vmax.*c);
Vmax_ac_denominator = (Km+c);
Vmax_ac = Vmax_ac_Numerator./Vmax_ac_denominator;

Fluxace(ace_points,1) = 0;
mu(ace_points,1) = 0;
Fluxfe3(ace_points,1) = 0;
Fluxhco3(ace_points,1) = 0;
Fluxh(ace_points,1) = 0;
% Fluxall = cell(2, ace_points);
   
for i=1:ace_points 
   % restoreEnvironment(environment);   
       p=strjoin(string([x(i) y(i)]));
    model = changeRxnBounds(model,'EX_ac(e)',[-Vmax_ac(i)],'l');
    sol=solveCobraLP(model);
 
    if sol.stat ~= 1 
                Position_new{1, i} = p;
                Fluxall_new{i} = zeros(fluxsize); 
    else       
                Fluxace(i)=sol.full(id_ace);
                Fluxfe3(i)=sol.full(id_fe3);
                Fluxhco3(i)=sol.full(id_hco3);
                Fluxh(i)=sol.full(id_h);
                mu(i)= sol.obj;
                Position_new{1, i} = p;
                Fluxall_new{i} = sol.full;
     end
    
end

%% to save time step results in global variables for further processing

if n==10
    Fluxall{n}=[Fluxall_new];
end

if n==20
    Fluxall{n}=[Fluxall_new];
end

if n==30
    Fluxall{n}=[Fluxall_new];
end
if n==40
    Fluxall{n}=[Fluxall_new];
end

if n==50
    Fluxall{n}=[Fluxall_new];
end

if n==80
    Fluxall{n}=[Fluxall_new];
end

if n==100
    Fluxall{n}=[Fluxall_new];
end


Vace=Fluxace;
Vhco3=Fluxhco3;
Vfe3=Fluxfe3;
Vh=Fluxh;
growth=mu;

   

% Vace_xy=[Vace cv_compartments1 x y t];

% 1 mmol/g/h=1mol/(1000g×3600s)×(96485C/mol)=0.0268A/g
% (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3673867/)

% GS_Currents_Fe3=-(0.0268)*(1e+06)*(8*Vace.*c2)*(VOLUME./ace_points);
% GS_Currents_Fe3_positions=[GS_Currents_Fe3 x y t];
% GS_Currents_Fe3_time=[sum(GS_Currents_Fe3) t(1)];

% fname = sprintf('C:\\Users\\acade\\Desktop\\Biofilm_RT\\v013\\MmodelSimDataPhase2\\M_model%s.mat', num2str(t(1)));   
% save(fname);

Acetate_phase2(:,n)=c;
GS_Currents_Fe3_phase2(:,n)=(0.0268)*(1e+06)*(Vfe3.*c2)*(VOLUME./ace_points);
% GS_Currents_Fe3(:,n)=-(0.0268)*(1e+06)*(8*Vace.*c2)*(VOLUME./ace_points);
% GS_Currents_Fe3_positions=[GS_Currents_Fe3 x y t];
GS_Currents_Fe3_phase2_time_series(n)=sum(GS_Currents_Fe3_phase2(:,n));
Acetate_phase2_time_series(n)=sum(Acetate_phase2(:,n))/ace_points;

end