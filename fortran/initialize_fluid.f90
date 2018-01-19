subroutine initialize_fluid
    use all_parameters
    implicit none
         
    integer :: i,j,k,counter_fl,check_fl
    double precision  :: ran1,x12,y12,z12,ke1_fluid,zz1, r
    double precision  :: avr_vel_fl_x,avr_vel_fl_y,avr_vel_fl_z,tx,ty,tz
    double precision  :: average_vel_fl_x,average_vel_fl_y,average_vel_fl_z
          
    counter_fl = 0
    avr_vel_fl_x = 0.0d0; avr_vel_fl_y = 0.0d0; avr_vel_fl_z = 0.0d0
          
    do while(counter_fl < no_of_fluid)
        tx = ran1(zzzz)*llx
        ty = ran1(zzzz)*lly
        tz = ran1(zzzz)*llz
        check_fl =1
        
        do j=1,no_of_colloid
            x12 = tx - pos_colloid(3*j-2)
            y12 = ty - pos_colloid(3*j-1)
            z12 = tz - pos_colloid(3*j)
            x12 = x12 - llx*anint(x12*inv_llx)
            y12 = y12 - lly*anint(y12*inv_lly)
            z12 = z12 - llz*anint(z12*inv_llz)

            r = dsqrt(x12**2+y12**2+z12**2)
            if(r.lt.sigma*0.5) then                 !fluid particle is at least one colloid_radius apart
                check_fl =0                             !from the colloid's centre
            endif
        enddo
        if(check_fl==1) then
            counter_fl = counter_fl +1
            pos_fl(3*counter_fl-2) = tx          !initialize fluid position
            pos_fl(3*counter_fl-1) = ty
            pos_fl(3*counter_fl) = tz
        endif
    enddo

    do j=1,no_of_fluid
        vel_fl(3*j-2) = (ran1(zzzz)-0.5d0)*vscale_fluid
        vel_fl(3*j-1) = (ran1(zzzz)-0.5d0)*vscale_fluid
        vel_fl(3*j)   = (ran1(zzzz)-0.5d0)*vscale_fluid
        avr_vel_fl_x  = avr_vel_fl_x + vel_fl(3*j-2)
        avr_vel_fl_y = avr_vel_fl_y + vel_fl(3*j-1)
        avr_vel_fl_z = avr_vel_fl_z + vel_fl(3*j)
    enddo

    avr_vel_fl_x = avr_vel_fl_x/dfloat(no_of_fluid)
    avr_vel_fl_y = avr_vel_fl_y/dfloat(no_of_fluid)
    avr_vel_fl_z = avr_vel_fl_z/dfloat(no_of_fluid)

    average_vel_fl_x = 0.0d0
    average_vel_fl_y = 0.0d0
    average_vel_fl_z = 0.0d0

    ke1_fluid=0.0d0
    do j=1,no_of_fluid
        vel_fl(3*j-2) = vel_fl(3*j-2) - avr_vel_fl_x
        vel_fl(3*j-1) = vel_fl(3*j-1) - avr_vel_fl_y
        vel_fl(3*j)   = vel_fl(3*j) - avr_vel_fl_z
        average_vel_fl_x = average_vel_fl_x + vel_fl(3*j-2)
        average_vel_fl_y = average_vel_fl_y + vel_fl(3*j-1)
        average_vel_fl_z = average_vel_fl_z + vel_fl(3*j)
        ke1_fluid = ke1_fluid + 0.50d0*mass_fl*(vel_fl(3*j-2)**2 + vel_fl(3*j-1)**2 + vel_fl(3*j)**2)
    enddo
end subroutine initialize_fluid
