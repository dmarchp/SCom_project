echo Generando script gnuplot para $1
sed "s/template/$1/" figures_template.gnu > figures_$1.gnu
