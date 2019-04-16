% Design of Square Footing

clc;
clear all;
format short g;

load value.mat
load input.mat
disp ("Design of Square Footing")
disp("\n")
%Compute the Area of the Footing
Load_on_Column = Load;
printf("Load_on_Column = %d N \n",Load_on_Column)

Assume_the_Self_Weight_of_Footing = (Self_Weight);
printf("Assume_the_Self_Weight_of_Footing = %d %% \n",Assume_the_Self_Weight_of_Footing)

Load_inc = (Self_Weight/100)*Load;
Self_Weight_of_Footing = Load_inc;
printf("Self_Weight_of_Footing = %d N \n",Self_Weight_of_Footing)

Vertical_Load_on_Column = Load+Load_inc;
printf("Vertical_Load_on_Column = %d N \n",Vertical_Load_on_Column)

disp("\n")
Required_Area_of_Footing = (Vertical_Load_on_Column/Soil_Pressure);
printf("Required_Area_of_Footing = %d mm^2 \n",Required_Area_of_Footing)

% Side of Rectangular Footing
Side_of_Footing = round(sqrt(Required_Area_of_Footing)/100)*100;
printf("Side_of_Footing = %d mm \n",Side_of_Footing)

Provided_Area_of_Footing = Side_of_Footing*Side_of_Footing;
printf("Provided_Area_of_Footing = %d mm^2 \n",Provided_Area_of_Footing)

% Net Upward pressure
Net_Upward_Pressure = (1.5*Load)/(Side_of_Footing*Side_of_Footing);
printf("Net_Upward_Pressure = %d N/mm^2 \n",Net_Upward_Pressure)

% Depth on the basis of moment
disp("\n")
disp("Depth on the basis of Bending Compression")


M = Net_Upward_Pressure*(Side_of_Footing/8)*(Side_of_Footing-Side_of_Column)*(Side_of_Footing-Side_of_Column);
printf("Moment_at_the_Face_of_the_Column = %d Nmm \n",M)

Factored_Moment_at_the_Face_of_the_Column = 1.5*M;
printf("Factored_Moment_at_the_Face_of_the_Column = %d Nmm \n",Factored_Moment_at_the_Face_of_the_Column)

xu_max_by_d = (700)/(1100+0.87*Fy);
Ru = (0.36*Fck*xu_max_by_d*(1-0.416*xu_max_by_d));

Effective_Depth_on_the_basis_of_Moment= round(sqrt((Factored_Moment_at_the_Face_of_the_Column)/(Ru*Side_of_Footing))/10)*10;
printf("Effective_Depth_on_the_basis_of_Moment = %d mm \n",Effective_Depth_on_the_basis_of_Moment)

Overall_Depth_on_the_basis_of_Moment = Effective_Depth_on_the_basis_of_Moment+Clear_Cover;
printf("Overall_Depth_on_the_basis_of_Moment  = %d mm \n",Overall_Depth_on_the_basis_of_Moment)

disp("\n")
disp("Depth on the basis of One Way Shear:-")
if(Fck>40)
fck=40;
elseif(Fck<40)
fck=Fck;
endif
Permissible_Shear_Stress1 = interp2(tables,tables,tables,fck,pt);
if(Overall_Depth_on_the_basis_of_Moment>300)
Overall_Depth_on_the_basis_of_Moment=300;
endif
k = interp1(table(:,1),table(:,2),Overall_Depth_on_the_basis_of_Moment);

Permissible_Shear_Stress = k*Permissible_Shear_Stress1;
printf("Permissible_Shear_Stress = %d N/mm^2 \n",Permissible_Shear_Stress) 

L_bey_col = (Side_of_Footing-Side_of_Column)/2;
d_one_way_shear = round((Net_Upward_Pressure*L_bey_col)/(Net_Upward_Pressure+Permissible_Shear_Stress));
printf("Eff_depth_one_way_shear = %d mm \n",d_one_way_shear)

disp("\n")
disp("Depth on the basis of Two Way Shear:-")
Width = Side_of_Column+d_one_way_shear;
Shear_Force_around_the_section = (Net_Upward_Pressure*(Side_of_Footing*Side_of_Footing-(Width*Width)));
printf("Shear_Force_around_the_section = %d N \n",Shear_Force_around_the_section)

% beta_c is the ratio of short side to lonf side of the column so for square column beta_c = 1
beta_c = 1;
Ks = (0.5+beta_c);
if(Ks>1)
Ks=1;
elseif(Ks<1)
Ks=Ks;
endif
tc = 0.25*sqrt(Fck);
depth_two_way_shear = [4*Ks*tc 4*Ks*tc*Side_of_Column -Shear_Force_around_the_section];
d_two_way_shear = max(roots(depth_two_way_shear));
printf("Eff_depth_two_way_shear = %d mm \n",d_two_way_shear)

disp("\n")
depth = [Effective_Depth_on_the_basis_of_Moment d_one_way_shear d_two_way_shear];
Eff_depth_of_footing = max(depth);
printf("Eff_depth_of_footing = %d mm \n",Eff_depth_of_footing)
Overall_depth_of_footing = round((Eff_depth_of_footing+Clear_Cover+dia)/10)*10;
printf("Overall_depth_of_footing = %d mm \n",Overall_depth_of_footing)

disp("\n")
% Design of Reinforcement
disp("Design of Steel Reinforcement:-")
Area_of_steel = 0.5*(Fck/Fy)*(1-sqrt(1-((4.6*M)/(Fck*Side_of_Footing*Eff_depth_of_footing*Eff_depth_of_footing))))*Side_of_Footing*Eff_depth_of_footing;
printf("Area_of_steel = %d mm^2 \n", Area_of_steel)

Area_of_one_bar = (pi/4)*(dia*dia);                                             
printf("Area_of_one = %d mm^2 \n", Area_of_one_bar)

No_of_Reinforced_Bars = ceil(Area_of_steel/Area_of_one_bar)
disp("Hence provide these number of Reinforcement Bars uniformuly distributed in each direction.")

disp("\n")
disp("Check for Development Length:-")
tbd = interp1 (mat(:,1),mat(:,2),Fck);
Development_Length = round((0.87*Fy)/(4*tbd*1.6))*dia;
printf("Development_Length = %d mm \n", Development_Length)
Length_Available_of_Bars = round((0.5*(Side_of_Footing - Side_of_Column) - (Side_Cover)));
printf("Length_Available_of_Bars = %d mm \n", Length_Available_of_Bars)

if(Development_Length<Length_Available_of_Bars)
disp("Hence Safe in Development Length")
elseif(Development_Length>Length_Available_of_Bars)
disp("Hence not Safe in Development Length")
endif

disp("\n")
% Check for Transfer of Load at the base
disp("Check for Transfer of Load at the base:-")
A2 = (Side_of_Column)*(Side_of_Column);

A1 = ((Side_of_Column)+2*(2*Overall_depth_of_footing))*((Side_of_Column)+2*(2*Overall_depth_of_footing));

square_root_A1_by_A2 = sqrt((A1)/(A2));
if(square_root_A1_by_A2>2)
square_root_A1_by_A2 = 2;
elseif(square_root_A1_by_A2<2)
square_root_A1_by_A2 = square_root_A1_by_A2;
endif

Permissible_Bearing_Stress = 0.45*Fck*square_root_A1_by_A2;
printf("Permissible_Bearing_Stress = %d N/mm^2 \n",Permissible_Bearing_Stress)
Actual_Bearing_Stress = (Load_on_Column)/(Side_of_Column*Side_of_Column);
printf("Actual_Bearing_Stress = %d N/mm^2 \n",Actual_Bearing_Stress)

if(Permissible_Bearing_Stress>Actual_Bearing_Stress)
disp("Hence Safe in Bearing Stress. Thus no Seperate Dowel Bars are required for the transfer of load. However, it is advisable to continue all the bars of the Column, into the Foundation.")
elseif(Permissible_Bearing_Stress<Actual_Bearing_Stress)
disp("Hence not Safe in Bearing Stress")
endif