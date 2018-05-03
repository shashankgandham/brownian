for colloid in {1..10}
do
	sed -i '2s/.*/EXECUTABLE\=a'$colloid'.out/' Makefile
    	sed -i '27s/.*/integer,parameter::no_of_colloid\='$colloid'/' define_parameter.f90
	make clean
	make
done
