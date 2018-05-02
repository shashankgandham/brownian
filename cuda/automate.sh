for density in 5 10 20 
do
	for colloid in {1..10} 
	do
		echo "nvprof --log-file vprof$density$colloid.txt ./a.out $density $colloid"
		nvprof --log-file vprof$density$colloid.txt ./a.out $density $colloid 
		rm -rf /tmp/.nvprof/*
	done
done
