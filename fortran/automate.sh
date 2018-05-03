sleep 300
for density in 5
do
	for colloid in {1..10} 
	do
		echo "Running $density $colloid"
		sudo nice -n -20 ./a$colloid.out > prof$density$colloid.txt
		gprof gmon.out > gmon$density$colloid.txt
		cat gmon$density$colloid.txt >> prof$density$colloid.txt
	done
done
