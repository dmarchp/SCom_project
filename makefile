objects=main.o read_net.o net_props.o sir.o
dep_objects=read_net.o net_props.o sir.o
mods=set_up_net.mod net_props.mod SIR.mod
compiler=gfortran-9
opt = -fbounds-check

main.exe : $(objects)
	mkdir -p results
	$(compiler) -o main.exe $(opt) $(objects)

$(mods) : $(dep_objects)

read_net.o : read_net.f90
	$(compiler) -c $(opt) read_net.f90

net_props.o: net_props.f90
	$(compiler) -c $(opt) net_props.f90
	
sir.o : sir.f90 read_net.f90
	$(compiler) -c $(opt) sir.f90

main.o : main.f90 $(mods)
	$(compiler) -c $(opt) main.f90

.PHONY: clean
clean:
	rm -f *.o
	rm -f *.mod
	
.PHONY: results
results:
	make
	./main.exe net$(net_size).dat ${seed}
	mkdir -p results/net$(net_size)/
	mv max_inf_infrec_lambda.dat results/net$(net_size)/
	mv pop_frac_evo_lambda_*.dat results/net$(net_size)/

