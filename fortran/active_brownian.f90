       program colloid
       use all_parameters
       implicit none
       real*8 :: zz1,ran1,random_expo,energy_colloid,ke_colloid2
       real*8 ::distance_x,distance_y,distance_z,distance
       real*8 ::mom1_x,mom1_y,mom1_z,mom2_x,mom2_y,mom2_z
       integer::j,i,l,mm,nnn
       character(len=30) :: charac_a,charac_b,charac_c,charac_d
       
       charac_b = 'snapshot_fluid'
!       charac_d = 'snapshot_colloid'
       charac_d = 'snapshot'
       
       open(39,file='cluster.dat')
       open(1111,file='final_colloid_position.dat')
       open(1112,file='final_colloid_velocity.dat')
       open(1113,file='final_fluid_position.dat')
       open(1114,file='final_fluid_velocity.dat')
       open(1116,file='ke_fluid.dat')
       open(1117,file='ke_colloid.dat')
        zz1 = ran1(zzz)
      if (file_colloid.eq.1) then !read from file
        
        do mm=1,no_of_colloid
         read(1111,fmt='(9g25.15)') pos_colloid(3*mm-2),pos_colloid(3*mm-1),pos_colloid(3*mm)
         read(1112,fmt='(9g25.15)') vel_colloid(3*mm-2),vel_colloid(3*mm-1),vel_colloid(3*mm)
        enddo
        do  mm=1,no_of_fluid
            read(1113,*)pos_fl(3*mm-2),pos_fl(3*mm-1),pos_fl(3*mm)
            read(1114,*)vel_fl(3*mm-2),vel_fl(3*mm-1),vel_fl(3*mm)
        enddo
        else
     call initialize_colloid
!       call initialize_activity
     call initialize_fluid
       endif
     !           write(*,*) 'INITIAL KE1',ke_colloid1*mass_colloid/(2.0d0*dfloat(no_of_colloid))
     !           write(*,*) 'INITIAL ANGULAR KE',ang_ke_colloid1*I_colloid/(2.0d0*dfloat(no_of_colloid))
     call create_box
!     write(*,*) maxval(box_neigh)
     call neighbour_list_mpcd
     call neighbour_list_md
     call compute_force_md
     call tumble
     write(*,*)"after tumble"   
!vel_fl=0.0d0
        do nn=1,niter
          write(*,*) nn
           
           call rotation_mpcd
                call run
           do l=1,n
              qq=(nn-1)*10+l
                    if(nn>90.and.nn>10000)       write(121,*) vel_colloid(1),vel_colloid(2),vel_colloid(3)
              old_force=f
!write(*,*) qq    
          call update_pos_md
!              old_ra=ra
               call neighbour_list_md
              call update_pos_mpcd    !update fluid position with MD time step for all fluid particles
               call neighbour_list_mpcd
              if(mod(qq,10).eq.0.and.nn>10000) call updown_velocity
              call fluid_colloid_collision
              call update_activity_direction
              call compute_force_md
              call update_velocity_colloid

              
           enddo
!           if(nn.gt.80000)then
!             call cluster
!             do i=1,n_cluster
!             write(39,*)nn,i,size_cluster(i)
!             enddo
!           endif

                ke_colloid=0.0d0;ke_fluid=0.0d0;ang_ke_colloid=0.0d0
           do   i=1,no_of_colloid
                ke_colloid = ke_colloid + vel_colloid(3*i-2)**2 + vel_colloid(3*i-1)**2 + vel_colloid(3*i)**2
           enddo
                ke_colloid = 0.50d0*mass_colloid*ke_colloid

           do   i=1,no_of_colloid
                ang_ke_colloid = ang_ke_colloid + ang_vel_colloid(3*i-2)**2 + ang_vel_colloid(3*i-1)**2 + ang_vel_colloid(3*i)**2
           enddo
                ang_ke_colloid=0.50d0*ang_ke_colloid*I_colloid
!              if(mod(nn,step_thermo)==0)  call thermostat_colloid

            write(1117,*) nn,ke_colloid,ang_ke_colloid   
                energy_colloid=potential_colloid+ke_colloid+ang_ke_colloid


           do   i=1,no_of_fluid
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

        write(115,*) mom_x,mom_y,mom_z

!       if(ke_colloid/dfloat(no_of_colloid).gt.30.0d0) then
!         write(8777,fmt='(9g25.15)') nn,pos_colloid(1),pos_colloid(2),pos_colloid(3)
         write(8778,fmt='(9g25.15)') nn,vel_colloid(1),vel_colloid(2),vel_colloid(3)
!         write(7778,*) pos_colloid(4),pos_colloid(5),pos_colloid(6)
!       endif
!         write(*,*)  nn, (energy_colloid+ke_fluid)  !total energy of colloid(except angular part)
!         write(5105,fmt='(9g25.15)')dfloat(nn),ke_fluid/dfloat(no_of_fluid), & 
!       ke_colloid/dfloat(no_of_colloid),ang_ke_colloid/dfloat(no_of_colloid)
!         write(94,*) vel_fl
!         write(034,*)vel_colloid
!      write(*,*) ke_colloid,nn

!      charac_a = charac_b
!      call addnumtostring(charac_a,nn)

!      open(071,file=charac_a,status='unknown')

!      do nnn=1,no_of_fluid
!        write(71,*) pos_fl(3*nnn-2),pos_fl(3*nnn-1),pos_fl(3*nnn)
!      end do

!      close(71)
!      charac_a = charac_b    !for fluid
      charac_c = charac_d
!      if(mod(nn,10).eq.0) then
!      call addnumtostring(charac_c,nn)
!      call addnumtostring(charac_a,nn) !for fluid
!      open(082,file=charac_a,status='unknown')  !for fluid
!      write(82,fmt='(9g25.15)') 0.5*lx,0.5*ly,0.5*lz,mom1_x,mom1_y,mom1_z  !for fluid
!      close(82)

!      open(072,file=charac_c,status='unknown')

!      do nnn=1,no_of_colloid
 
!       write(72,fmt='(9g25.15)') pos_colloid(3*nnn-2),pos_colloid(3*nnn-1),pos_colloid(3*nnn), &
!     & vel_colloid(3*nnn-2),vel_colloid(3*nnn-1),vel_colloid(3*nnn)
!      end do

!      close(72)
!      endif

        enddo !niter


!        do  mm=1,no_of_fluid
!            write(1113,*)pos_fl(3*mm-2),pos_fl(3*mm-1),pos_fl(3*mm)
!            write(1114,*)vel_fl(3*mm-2),vel_fl(3*mm-1),vel_fl(3*mm)
!        enddo
        do mm=1,no_of_colloid
         write(1111,fmt='(9g25.15)') pos_colloid(3*mm-2),pos_colloid(3*mm-1),pos_colloid(3*mm)
         write(1112,fmt='(9g25.15)') vel_colloid(3*mm-2),vel_colloid(3*mm-1),vel_colloid(3*mm)
        enddo
      end program colloid
