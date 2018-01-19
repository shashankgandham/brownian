            subroutine initialize_activity
            use all_parameters
            implicit none
            integer::i
            double precision :: aaa,ran1
			double precision x, y, z
            !!!!!!!FIND THE DIRECTION OF THE ACTIVE VELOCITY
            do i=1,no_of_colloid
              ra(3*i-2)= pos_colloid(3*i-2)-ran1(zzzz)*llx
              ra(3*i-1)= pos_colloid(3*i-1)-ran1(zzzz)*lly
              ra(3*i)= pos_colloid(3*i)-ran1(zzzz)*llz
              x = ra(3*i - 2)
			  y = ra(3*i - 1)
			  z = ra(3*i)
              aaa=dsqrt(x**2+y**2+z**2)

              ra(3*i-2)=ra(3*i-2)/aaa !RA IS THE UNIT VECTOR TOWARDS THE ACTIVE VELOCITY
              ra(3*i-1)=ra(3*i-1)/aaa
              ra(3*i)=ra(3*i)/aaa

              vel_colloid(3*i-2)=vel_colloid(3*i-2)+ra(3*i-2)*v0
              vel_colloid(3*i-1)=vel_colloid(3*i-1)+ra(3*i-1)*v0
              vel_colloid(3*i)=vel_colloid(3*i)+ra(3*i)*v0

            enddo
            end subroutine initialize_activity



