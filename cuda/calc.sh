calculate(){
	temp=$1
	temp=${temp::-1}
	if [[ "$temp" == *m ]]; then 
		temp=${temp::-1}
		temp=$(echo "scale=8; $temp/1000" | bc)
	fi
	if [[ "$temp" == *u ]]; then 
		temp=${temp::-1}
		temp=$(echo "scale=8; $temp/1000000" | bc)
	fi
	echo $temp
}

echo > timings.txt
for density in 5 10 20
do
	for colloid in 1 10
	do
		echo "" >> timings.txt
		time=0
		echo "$density $colloid" >> timings.txt
		output=$(grep neighbour_list_mpcd prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		output=$(grep sieve prof$density$colloid.txt)
		time2=$(echo $output | awk '{print $4}')
		time2=$(calculate $time2)
		output=$(grep boxpart prof$density$colloid.txt)
		time3=$(echo $output | awk '{print $2}')
		time3=$(calculate $time3)
		time=$(echo "$time1 + $time2 + $time3" | bc)
		echo "neighbour_list_mpcd $time" >> timings.txt

		output=$(grep cellpart prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		output=$(grep cellvel prof$density$colloid.txt)
		time2=$(echo $output | awk '{print $4}')
		time2=$(calculate $time2)
		output=$(grep velfl prof$density$colloid.txt)
		time3=$(echo $output | awk '{print $2}')
		time3=$(calculate $time3)
		output=$(grep rotate_mat prof$density$colloid.txt)
		time4=$(echo $output | awk '{print $4}')
		time4=$(calculate $time4)
		output=$(grep rotate\( prof$density$colloid.txt)
		time5=$(echo $output | awk '{print $2}')
		time5=$(calculate $time5)
		time=$(echo "$time1 + $time2 + $time3 + $time4 + $time5" | bc)
		echo "rotation_mpcd $time" >> timings.txt


		time=0
		output=$(grep fluid_colloid_collision prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		output=$(grep reduce prof$density$colloid.txt)
		time2=$(echo $output | awk '{print $4}')
		time2=$(calculate $time2)
		output=$(grep d_dump prof$density$colloid.txt)
		time3=$(echo $output | awk '{print $2}')
		time3=$(calculate $time3)
		output=$(grep update_fcc prof$density$colloid.txt)
		time4=$(echo $output | awk '{print $2}')
		time4=$(calculate $time4)
		time=$(echo "$time1 + $time2 + $time3 + $time4" | bc)
		echo "fluid_colloid_collision $time" >> timings.txt
		
		time=0
		output=$(grep update_pos_mpcd prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		time=$(echo "$time1" | bc)
		echo "update_pos_mpcd $time" >> timings.txt

		time=0
		output=$(grep velc prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		output=$(grep nbrc prof$density$colloid.txt)
		time2=$(echo $output | awk '{print $4}')
		time2=$(calculate $time2)
		time=$(echo "$time1 + $time2" | bc)
		echo "run $time" >> timings.txt

		time=0
		output=$(grep helper_upd prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		output=$(grep calc_upd prof$density$colloid.txt)
		time2=$(echo $output | awk '{print $4}')
		time2=$(calculate $time2)
		time=$(echo "$time1 + $time2" | bc)
		echo "updown_velocity $time" >> timings.txt
	
		time=0
		output=$(grep update_vel_colloid prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		time=$(echo "$time1" | bc)
		echo "update_velocity_colloid $time" >> timings.txt

		time=0
		output=$(grep create_box prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		time=$(echo "$time1" | bc)
		echo "create_box $time" >> timings.txt

		time=0
		output=$(grep compute_force_md prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		time=$(echo "$time1" | bc)
		echo "compute_force_md $time" >> timings.txt

		time=0
		output=$(grep neighbour_list_md prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		time=$(echo "$time1" | bc)
		echo "neighbour_list_md $time" >> timings.txt

		time=0
		output=$(grep initialize_fluid prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		time=$(echo "$time1" | bc)
		echo "initialize_fluid $time" >> timings.txt

		time=0
		output=$(grep update_activity_direction prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		time=$(echo "$time1" | bc)
		echo "update_activity_direction $time" >> timings.txt

		time=0
		output=$(grep tumble prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		time=$(echo "$time1" | bc)
		echo "tumble $time" >> timings.txt

		time=0
		output=$(grep update_pos_md prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		time=$(echo "$time1" | bc)
		echo "update_pos_md $time" >> timings.txt

		time=0
		output=$(grep initialize_colloid prof$density$colloid.txt)
		time1=$(echo $output | awk '{print $2}')
		time1=$(calculate $time1)
		time=$(echo "$time1" | bc)
		echo "initialize_colloid $time" >> timings.txt



	done
done
