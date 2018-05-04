for density in 5 10 20 
do
    	sed -i '22s/.*/integer,parameter::no_of_fluid\=lx*ly*lz*'$density'/' define_parameter.f90
	for colloid in {1..10}
	do
		sed -i '2s/.*/EXECUTABLE\=a'$colloid'.out/' Makefile
    		sed -i '27s/.*/integer,parameter::no_of_colloid\='$colloid'/' define_parameter.f90
		make clean
		make
	done
done
