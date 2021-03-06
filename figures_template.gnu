
set term pngcairo enhanced
set output "infected_time_evo_template.png"

set ylabel "f_i"
set xlabel "t / {/Symbol d}"
set key title "{/Symbol l}"
set yrange[0:1]
set xrange[0:6]

file1 = "results/nettemplate/pop_frac_evo_lambda_1.dat"
file2 = "results/nettemplate/pop_frac_evo_lambda_3.dat"
file3 = "results/nettemplate/pop_frac_evo_lambda_5.dat"
file4 = "results/nettemplate/pop_frac_evo_lambda_8.dat"
file5 = "results/nettemplate/pop_frac_evo_lambda_10.dat"
file6 = "results/nettemplate/pop_frac_evo_lambda_20.dat"
file7 = "results/nettemplate/pop_frac_evo_lambda_30.dat"
file8 = "results/nettemplate/pop_frac_evo_lambda_40.dat"
file9 = "results/nettemplate/pop_frac_evo_lambda_50.dat"
file10 = "results/nettemplate/pop_frac_evo_lambda_60.dat"
file11 = "results/nettemplate/pop_frac_evo_lambda_70.dat"
file12 = "results/nettemplate/pop_frac_evo_lambda_80.dat"
file13 = "results/nettemplate/pop_frac_evo_lambda_90.dat"
file14 = "results/nettemplate/pop_frac_evo_lambda_100.dat"

plot file1 i 1 u 2:4 w l lc rgb "purple" t"0.1", \
file2 i 0 u 2:4 w l lc rgb "dark-blue" t"0.3", \
file3 i 0 u 2:4 w l lc rgb "web-blue" t"0.5", \
file4 i 0 u 2:4 w l lc rgb "cyan" t"0.8", \
file5 i 0 u 2:4 w l lc rgb "green" t"1.0", \
file6 i 0 u 2:4 w l lc rgb "greenyellow" t"2.0", \
file7 i 0 u 2:4 w l lc rgb "goldenrod" t"3.0", \
file8 i 0 u 2:4 w l lc rgb "red" t"4.0", \
file9 i 0 u 2:4 w l lc rgb "dark-red" t"5.0";
#file10 i 0 u 2:4 w l lc rgb "pink" t"6.0", \
#file11 i 0 u 2:4 w l lc rgb "sienna1" t"7.0", \
#file12 i 0 u 2:4 w l lc rgb "bisque" t"8.0", \
#file13 i 0 u 2:4 w l lc rgb "gray0" t"9.0", \
#file14 i 0 u 2:4 w l lc rgb "gray70" t"10.0";


set output "max_infected_lambda_template.png"
set ylabel "max f_i"
set xlabel "{/Symbol l} / {/Symbol d}"
unset key

plot "results/nettemplate/max_inf_infrec_lambda.dat" i 0 u 1:2 w l notitle;
