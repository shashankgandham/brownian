program colloid
    use all_parameters
	implicit none
	double precision ::mom1_x,mom1_y,mom1_z,mom2_x,mom2_y,mom2_z,energy_colloid
	integer::j,i,l

	call initialize_colloid
    call initialize_fluid
    call create_box
	call neighbour_list_mpcd
	call neighbour_list_md
	call compute_force_md
	call tumble
    write(*,*) "After Tumble"
    do nn = 1, niter
		write(*,*) nn
		call rotation_mpcd
		call run
        do l=1,n
			qq=(nn-1)*10+l
			old_force=f
			call update_pos_md
			call neighbour_list_md
			call update_pos_mpcd  
			call neighbour_list_mpcd
			if(mod(qq,10).eq.0.and.nn>10000) call updown_velocity
		    call fluid_colloid_collision
            call update_activity_direction
			call compute_force_md
			call update_velocity_colloid
		enddo
        ke_colloid=0.0d0;ke_fluid=0.0d0;ang_ke_colloid=0.0d0

		do   i=1,no_of_colloid
        	ke_colloid = ke_colloid + vel_colloid(3*i-2)**2 + vel_colloid(3*i-1)**2 + vel_colloid(3*i)**2
       	enddo
     	ke_colloid = 0.50d0*mass_colloid*ke_colloid

		do i=1,no_of_colloid
        	ang_ke_colloid = ang_ke_colloid + ang_vel_colloid(3*i-2)**2 + ang_vel_colloid(3*i-1)**2 + ang_vel_colloid(3*i)**2
		enddo
        ang_ke_colloid=0.50d0*ang_ke_colloid*I_colloid
		energy_colloid=potential_colloid+ke_colloid+ang_ke_colloid

		do i=1,no_of_fluid
			ke_fluid = ke_fluid + vel_fl(3*i-2)**2 + vel_fl(3*i-1)**2 + vel_fl(3*i)**2
        enddo
        ke_fluid=0.50d0*ke_fluid*mass_fl

		mom1_x=0.0d0;mom1_y=0.0d0;mom1_z=0.0d0
        mom2_x=0.0d0;mom2_y=0.0d0;mom2_z=0.0d0

        do i=1,no_of_fluid
			mom1_x=mom1_x+mass_fl*vel_fl(3*i-2)
			mom1_y=mom1_y+mass_fl*vel_fl(3*i-1)
			mom1_z=mom1_z+mass_fl*vel_fl(3*i)
		enddo

		do j=1,no_of_colloid
			mom2_x=mom2_x+mass_colloid*vel_colloid(3*j-2)
			mom2_y=mom2_y+mass_colloid*vel_colloid(3*j-1)
			mom2_z=mom2_z+mass_colloid*vel_colloid(3*j)
		enddo
		mom_x=mom1_x+mom2_x
		mom_y=mom1_y+mom2_y
		mom_z=mom1_z+mom2_z
	enddo
end program colloid
