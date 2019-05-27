# plot file for diff_eq_apr_25
#***************************************************************************
# FIRST PLOT
#***************************************************************************

set term x11 1 

#___________________________________________________________________________
# PLOT STYLE

set title\
"Small angle approximation"

# set xlabel 'x'
# set ylabel 'y'

# set key at 0.16,-3.3
# set key left

# set xrange[-3.5:-1]
# set yrange[-16:2]

# set timestamp

#___________________________________________________________________________
# FIT THE CURVES
# 
# f1(x) = a1*x + b1
# fit[-2.5:-1] f1(x) "diffeq_test2.dat" using (log10($1)) : (log10( abs($2 - $4) / $4)) via a1, b1
# fit1_title = sprintf("%-+4.1f*x %-+4.1f",a1,b1)
#
# f2(x) = a2*x + b2
# fit[-2.5:-1.5] f2(x) "diffeq_test2.dat" using (log10($1)) : (log10( abs($3 - $4) / $4)) via a2, b2
# fit2_title = sprintf("%-+4.1f*x %-+4.1f",a2,b2)
#
#___________________________________________________________________________
# PLOT

plot  "diff_eq_apr_27.dat" using 1:2 title 'theta(t)'\
     ,"diff_eq_apr_27.dat" using 1:3 title 'theta_dot(t)'\
#     ,a1*x + b1 title fit1_title\
#    ,"diffeq_test2.dat" using (log10($1)):(log10(abs($3 - $4) / $4)) title 'Runge-Kutta-4'\

	 
	 

#***************************************************************************
# SECOND PLOT
#***************************************************************************

set term x11 2 title 'Phase Space'

#___________________________________________________________________________
# PLOT STYLE

set title 'State Space'

# set xlabel 'x'
# set ylabel 'y'

# set key left

# set xrange[0.00011:0.5]
# set yrange[y1:y2]

# set timestamp

#___________________________________________________________________________
# FIT THE CURVES
# 
# f1(x) = a1*x + b1
# fit[-2.5:-1] f1(x) "diffeq_test2.dat" using (log10($1)) : (log10( abs($2 - $4) / $4)) via a1, b1
# fit1_title = sprintf("%-+4.1f*x %-+4.1f",a1,b1)
#
# f2(x) = a2*x + b2
# fit[-2.5:-1.5] f2(x) "diffeq_test2.dat" using (log10($1)) : (log10( abs($3 - $4) / $4)) via a2, b2
# fit2_title = sprintf("%-+4.1f*x %-+4.1f",a2,b2)
#
#___________________________________________________________________________
# PLOT

plot  "diff_eq_apr_27.dat" using 2:3
#     ,a1*x + b1 title fit1_title\
#    ,"diffeq_test2.dat" using (log10($1)):(log10(abs($3 - $4) / $4)) title 'Runge-Kutta-4'\




