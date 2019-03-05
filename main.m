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
printf("Load_on_Column = %d KN \n",Load_on_Column)

Assume_the_Self_Weight_of_Footing = (Self_Weight);
printf("Assume_the_Self_Weight_of_Footing = %d %% \n",Assume_the_Self_Weight_of_Footing)

Load_inc = (Self_Weight/100)*Load;
Self_Weight_of_Footing = Load_inc;
printf("Self_Weight_of_Footing = %d KN \n",Self_Weight_of_Footing)

Vertical_Load_on_Column = Load+Load_inc;
printf("Vertical_Load_on_Column = %d KN \n",Vertical_Load_on_Column)

disp("\n")
Required_Area_of_Footing = (Vertical_Load_on_Column/Soil_Pressure);
printf("Required_Area_of_Footing = %d m^2 \n",Required_Area_of_Footing)

Provided_Area_of_Footing = round(Vertical_Load_on_Column/Soil_Pressure);
printf("Provided_Area_of_Footing = %d m^2 \n",Provided_Area_of_Footing)

% Side of Rectangular Footing

Side_of_Footing = round(sqrt(Provided_Area_of_Footing)/0.1)*0.1;
printf("Side_of_Footing = %d m^2 \n",Side_of_Footing)

% Net Upward pressure
Net_Upward_Pressure = (Load)/(Side_of_Footing*Side_of_Footing);
printf("Net_Upward_Pressure = %d KN/m^2 \n",Net_Upward_Pressure)

% Depth on the basis of moment
disp("\n")
disp("Depth on the basis of Bending Compression")


M = Net_Upward_Pressure*(Side_of_Footing/8)*(Side_of_Footing-Side_of_Column)*(Side_of_Footing-Side_of_Column)*1000000;
printf("Moment at the Face of the Column = %d Nmm \n",M)

Factored_Moment_at_the_Face_of_the_Column = 1.5*M;
printf("Factored_Moment_at_the_Face_of_the_Column = %d Nmm \n",Factored_Moment_at_the_Face_of_the_Column)

xu_max_by_d = (700)/(1100+0.87*Fy);
Ru = (0.36*Fck*xu_max_by_d*(1-0.416*xu_max_by_d));

Effective_Depth_on_the_basis_of_Moment= round(sqrt((Factored_Moment_at_the_Face_of_the_Column)/(Ru*Side_of_Footing*1000))/10)*10;
printf("Effective_Depth_on_the_basis_of_Moment = %d mm \n",Effective_Depth_on_the_basis_of_Moment)

Overall_Depth_on_the_basis_of_Moment = Effective_Depth_on_the_basis_of_Moment+Clear_Cover;
printf("Overall_Depth_on_the_basis_of_Moment  = %d mm \n",Overall_Depth_on_the_basis_of_Moment)

disp("\n")
disp("Depth on the basis of One Way Shear")
Permissible_Shear_Stress1 = interp2(tables,tables,tables,Fck,pt);
k = interp1(table(:,1),table(:,2),Overall_Depth_on_the_basis_of_Moment);

Permissible_Shear_Stress = k*Permissible_Shear_Stress1;
printf("Permissible_Shear_Stress = %d N/mm^2 \n",Permissible_Shear_Stress) 

Effective_Depth_on_the_basis_of_One_Way_Shear = round((((0.75*Net_Upward_Pressure*Side_of_Footing-0.75*Net_Upward_Pressure*Side_of_Column))/(Permissible_Shear_Stress+0.0015*Net_Upward_Pressure))/10)*10;
printf("Effective_Depth_on_the_basis_of_One_Way_Shear = %d mm \n",Effective_Depth_on_the_basis_of_One_Way_Shear) 
Overall_Depth_One_Way_Shear = Effective_Depth_on_the_basis_of_One_Way_Shear+Clear_Cover;
printf("Overall_Depth_on_the_basis_of_One_Way_Shear = %d mm \n",Overall_Depth_One_Way_Shear)

disp("\n")
if(Effective_Depth_on_the_basis_of_Moment>Effective_Depth_on_the_basis_of_One_Way_Shear)
Effective_Depth_of_the_Footing = Effective_Depth_on_the_basis_of_Moment;

Overall_Depth_of_the_Footing = Effective_Depth_of_the_Footing+Clear_Cover_for_Footing;
printf("Overall_Depth_of_the_Footing = %d mm \n",Overall_Depth_of_the_Footing)

Effective_Depth_of_Footing_in_one_direction = Overall_Depth_of_the_Footing-Clear_Cover;
printf("Effective_Depth_of_Footing_in_one_direction = %d mm \n",Effective_Depth_of_Footing_in_one_direction)

Effective_Depth_of_Footing_in_other_direction = Effective_Depth_of_Footing_in_one_direction-dia;
printf("Effective_Depth_of_Footing_in_other_direction = %d mm \n",Effective_Depth_of_Footing_in_other_direction)

elseif(Effective_Depth_on_the_basis_of_Moment<Effective_Depth_on_the_basis_of_One_Way_Shear)
Effective_Depth_of_the_Footing = Effective_Depth_on_the_basis_of_One_Way_Shear;

Overall_Depth_of_the_Footing = Effective_Depth_of_the_Footing+Clear_Cover_for_Footing;
printf("Overall_Depth_of_the_Footing = %d mm \n",Overall_Depth_of_the_Footing)

Effective_Depth_of_Footing_in_one_direction = Overall_Depth_of_the_Footing-Clear_Cover;
printf("Effective_Depth_of_Footing_in_one_direction = %d mm \n",Effective_Depth_of_Footing_in_one_direction)

Effective_Depth_of_Footing_in_other_direction = Effective_Depth_of_Footing_in_one_direction-dia;
printf("Effective_Depth_of_Footing_in_other_direction = %d mm \n",Effective_Depth_of_Footing_in_other_direction)

endif

disp("\n")
% Check for Two Way Shear
disp("Check for Two Way Shear")
Width = Side_of_Column*1000+Effective_Depth_of_the_Footing;
Shear_Force_around_the_section = (Net_Upward_Pressure*(Side_of_Footing*Side_of_Footing-(Width/1000*Width/1000)));
printf("Shear_Force_around_the_section = %d KN \n",Shear_Force_around_the_section)

Factored_Shear_Force = 1.5*Shear_Force_around_the_section;
printf("Factored_Shear_Force = %d KN \n",Factored_Shear_Force)

Actual_Shear_Stress = (Factored_Shear_Force*1000)/(4*Width*Effective_Depth_of_the_Footing);
printf("Actual_Shear_Stress = %d N/mm^2 \n",Actual_Shear_Stress)

Permissible_Shear_Stress2 = 0.25*sqrt(Fck);
printf("Permissible_Shear_Stress = %d N/mm^2 \n",Permissible_Shear_Stress2)

if(Actual_Shear_Stress<Permissible_Shear_Stress2)
disp("Hence Safe in Two Way Shear")
elseif(Actual_Shear_Stress>Permissible_Shear_Stress2)
disp("Hence not Safe in Two Way Shear")
endif

disp("\n")
% Design of Reinforcement
disp("Design of Steel Reinforcement")
Area_of_steel = 0.5*(Fck/Fy)*(1-sqrt(1-((4.6*Factored_Moment_at_the_Face_of_the_Column)/(Fck*Side_of_Footing*1000*Effective_Depth_of_Footing_in_other_direction*Effective_Depth_of_Footing_in_other_direction))))*Side_of_Footing*1000*Effective_Depth_of_Footing_in_other_direction;
printf("Area_of_steel = %d mm^2 \n", Area_of_steel)

Area_of_one_bar = (pi/4)*(dia*dia);                                             
printf("Area_of_one = %d mm^2 \n", Area_of_one_bar)

No_of_Reinforced_Bars = ceil(Area_of_steel/Area_of_one_bar)
disp("Hence provide these number of Reinforcement Bars uniformuly distributed in each direction.")

disp("\n")

disp("Check for Development Length")
tbd = interp1 (mat(:,1),mat(:,2),Fck_for_Column);

Development_Length = ceil((Fy)/(4*tbd*1.6))*dia;
printf("Development_Length = %d mm \n", Development_Length)
Length_Available_of_Bars = round((0.5*(Side_of_Footing*1000 - Side_of_Column*1000) - (Side_Cover)));
printf("Length_Available_of_Bars = %d mm \n", Length_Available_of_Bars)

if(Development_Length<Length_Available_of_Bars)
disp("Hence Safe in Development Length")
elseif(Development_Length>Length_Available_of_Bars)
disp("Hence not Safe in Development Length")
endif

disp("\n")
% Check for Transfer of Load at the base
disp("Check for Transfer of Load at the base")
A2 = (Side_of_Column*1000)*(Side_of_Column*1000);

A1 = ((Side_of_Column*1000)+2*(2*Overall_Depth_of_the_Footing))*((Side_of_Column*1000)+2*(2*Overall_Depth_of_the_Footing));


square_root_A1_by_A2 = sqrt((A1)/(A2));
if(square_root_A1_by_A2>2)
square_root_A1_by_A2 = 2;
elseif(square_root_A1_by_A2<2)
square_root_A1_by_A2 = square_root_A1_by_A2;
endif

Permissible_Bearing_Stress = 0.45*Fck*square_root_A1_by_A2;
printf("Permissible_Bearing_Stress = %d N/mm^2 \n",Permissible_Bearing_Stress)
Actual_Bearing_Stress = (1.5*Load_on_Column*1000)/(Side_of_Column*1000*Side_of_Column*1000);
printf("Actual_Bearing_Stress = %d N/mm^2 \n",Actual_Bearing_Stress)

if(Permissible_Bearing_Stress>Actual_Bearing_Stress)
disp("Hence Safe in Bearing Stress. Thus no Seperate Dowel Bars are required for the transfer of load. However, it is advisable to continue all the bars of the Column, into the Foundation.")
elseif(Permissible_Bearing_Stress<Actual_Bearing_Stress)
disp("Hence not Safe in Bearing Stress")
endif



