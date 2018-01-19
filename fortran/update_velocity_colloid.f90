          subroutine update_velocity_colloid
                 use all_parameters
               implicit none
               integer i
               double precision :: dtby2m

                dtby2m = dt/(mass_colloid*2.0d0)
!               mom_x=0.0d0; mom_y=0.0d0; mom_z=0.0d0 ! NOTE

              do i=1,no_of_colloid
                vel_colloid(3*i-2)=vel_colloid(3*i-2)+(old_force(3*i-2)+f(3*i-2))*dtby2m ! NOTE
                vel_colloid(3*i-1)=vel_colloid(3*i-1)+(old_force(3*i-1)+f(3*i-1))*dtby2m
                vel_colloid(3*i)=vel_colloid(3*i)+(old_force(3*i)+f(3*i))*dtby2m


!                mom_x = mom_x + mass_colloid*vel_colloid(3*i-2)
!                mom_y = mom_y + mass_colloid*vel_colloid(3*i-1)
!                mom_z = mom_z + mass_colloid*vel_colloid(3*i)
              enddo

!              write(*,*) mom_x,mom_y,mom_z,nn

!              write(*,*) ke_colloid


              end subroutine update_velocity_colloid
