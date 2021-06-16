objects=main.o read_net.o net_props.o
dep_objects=read_net.o net_props.o
mods=set_up_net.mod net_props.mod
compiler=gfortran-9
opt = -fbounds-check

main.x : $(objects)
	$(compiler) -o main.x $(opt) $(objects)

$(mods) : $(dep_objects)

read_net.o : read_net.f90
	$(compiler) -c $(opt) read_net.f90

net_props.o: net_props.f90
	$(compiler) -c $(opt) net_props.f90

main.o : main.f90 $(mods)
	$(compiler) -c $(opt) main.f90

.PHONY: clean
clean:
	rm -f *.o
	rm -f *.mod
