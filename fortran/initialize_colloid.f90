subroutine initialize_colloid
    use all_parameters
    implicit none
    integer:: i,j,k,counter,check,mm, nofp
    double precision:: ran1,x12,y12,z12,r,temp,tx,ty,tz
    double precision:: avr_vel_colloid_x,avr_vel_colloid_y,avr_vel_colloid_z
    double precision:: average_vel_colloid_x,average_vel_colloid_y,average_vel_colloid_z
      
    nofp=0
    do k=40,(lx-1)*10,50
        do j=40,(ly-1)*10,50
            do i=40,(lz-1)*10,50
                nofp=nofp+1
                if(nofp.le.no_of_colloid)then
                    pos_colloid(3*nofp-2)=dfloat(i)/10d0  
                    pos_colloid(3*nofp-1)=dfloat(j)/10d0  
                    pos_colloid(3*nofp)=dfloat(k)/10d0
                endif
            enddo
        enddo
    enddo

    counter = 0
    avr_vel_colloid_x = 0.0d0
    avr_vel_colloid_y = 0.0d0
    avr_vel_colloid_z = 0.0d0
    mm = 0    
    do while(counter < no_of_colloid)
        tx = ran1(zzzz)*llx
        ty = ran1(zzzz)*lly
        tz = ran1(zzzz)*llz
           mm = mm+1
        check =1
        do j=1,counter
            x12 = tx - pos_colloid(3*j-2)
            y12 = ty - pos_colloid(3*j-1)
            z12 = tz - pos_colloid(3*j)
            if(abs(x12)>llxby2) x12 = llx - abs(x12)
            if(abs(y12)>llyby2) y12 = lly - abs(y12)
            if(abs(z12)>llzby2) z12 = llz - abs(z12)
            r = dsqrt(x12**2+y12**2+z12**2)
            if(r.lt.space_limit) then
                check =0
            endif
        enddo
        if(check==1) then
            counter = counter+1
            pos_colloid(3*counter-2) = tx
            pos_colloid(3*counter-1) = ty
            pos_colloid(3*counter) = tz
        endif
    enddo
    do j=1,no_of_colloid
        tx = ran1(zzzz) - 0.5d0
        ty = ran1(zzzz) - 0.5d0
        tz = ran1(zzzz) - 0.5d0
        vel_colloid(3*j-2) =  tx*vscale_colloid
        vel_colloid(3*j-1) =  ty*vscale_colloid
        vel_colloid(3*j)   =  tz*vscale_colloid
        avr_vel_colloid_x = avr_vel_colloid_x + vel_colloid(3*j-2)
        avr_vel_colloid_y = avr_vel_colloid_y + vel_colloid(3*j-1)
        avr_vel_colloid_z = avr_vel_colloid_z + vel_colloid(3*j)
    enddo
    avr_vel_colloid_x = avr_vel_colloid_x/dfloat(no_of_colloid)
    avr_vel_colloid_y = avr_vel_colloid_y/dfloat(no_of_colloid)
    avr_vel_colloid_z = avr_vel_colloid_z/dfloat(no_of_colloid) 

    average_vel_colloid_x = 0.0d0
    average_vel_colloid_y = 0.0d0
    average_vel_colloid_z = 0.0d0

    ke_colloid1=0.0d0
    do j=1,no_of_colloid
        vel_colloid(3*j-2) = vel_colloid(3*j-2) - avr_vel_colloid_x
        vel_colloid(3*j-1) = vel_colloid(3*j-1) - avr_vel_colloid_y
        vel_colloid(3*j) = vel_colloid(3*j) - avr_vel_colloid_z
        average_vel_colloid_x = average_vel_colloid_x + vel_colloid(3*j-2)
        average_vel_colloid_y = average_vel_colloid_y + vel_colloid(3*j-1)
        average_vel_colloid_z = average_vel_colloid_z + vel_colloid(3*j)
        temp = vel_colloid(3*j-2)**2 + vel_colloid(3*j-1)**2 + vel_colloid(3*j)**2
        ke_colloid1 = ke_colloid1 + temp
    enddo
    do j=1,no_of_colloid
        ang_vel_colloid(3*j-2) = (ran1(zzzz)-0.5d0)*ang_vscale_colloid
        ang_vel_colloid(3*j-1) = (ran1(zzzz)-0.5d0)*ang_vscale_colloid
        ang_vel_colloid(3*j)   = (ran1(zzzz)-0.5d0)*ang_vscale_colloid
    enddo
    ang_ke_colloid1=0.0d0
    do j=1,no_of_colloid
        temp = ang_vel_colloid(3*j-2)**2 + ang_vel_colloid(3*j-1)**2 + ang_vel_colloid(3*j)**2
        ang_ke_colloid1 = ang_ke_colloid1 + temp
    enddo
end subroutine
