sleep 300
make clean && make
for density in 5
do
	for colloid in {1..10} 
	do
		echo "Running $density $colloid"
		sudo nice -n -20 ./a.out $density $colloid > prof$density$colloid.txt
		gprof a.out >> prof$density$colloid.txt
	done
done
