for density in 5 10 20 
do
	for colloid in {1..10} 
	do
		echo "Running $density $colloid"
		./a.out $density $colloid > cpprof$density$colloid.txt
		gprof gmon.out > gmon$density$colloid.txt
	done
done
