sleep 300
for density in 5 10 20 
do
	for colloid in {1..10} 
	do
		echo "nvprof --log-file vprof$density$colloid.txt ./a.out $density $colloid"
		sudo nice -n -20 ./a.out $density $colloid > prof$density$colloid.txt
		sudo nice -n -20 /usr/local/cuda-9.1/bin/nvprof --log-file vprof$density$colloid.txt ./a.out $density $colloid 
		cat vprof$density$colloid.txt >> prof$density$colloid.txt
		rm -rf /tmp/.nvprof/*
	done
done
