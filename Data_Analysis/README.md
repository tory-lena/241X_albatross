# 241X_albatross

################################
sdlog2_dump.py is from https://github.com/lafranchi/Firmware, which
converts log data to csv file
python sdlog2_dump.py log_name.px4log > log_output.csv

################################
Log_Analysis.py read csv file

1. Level flight
wind velcity = 0.0
Lift = weight
Thrust = Drag


The air velcity is  U = sqrt(GPOS_VelE^2 + GPOS_VelN^2 +GPOS_VelD^2)
angle of attack alpha = ATT_Pitch ??????
Power data: PWR_*
Lift = (1/2) rho * U^2 * S * CL(alpha)
Stall_Speed = sqrt(2*W/(S*rho*CL_max))
Plot alpha vs CL
Plot U vs diff PWR_*

 
2.0 Glide test
glide angle = alpha = arctan(VelD/sqrt(GPOS_VelE^2 + GPOS_VelN^2)) ????
D/L = CD/CL = tan(alpha)