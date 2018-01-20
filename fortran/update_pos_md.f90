subroutine update_pos_md
	use all_parameters
	implicit none
	integer i
	double precision :: ddt

	ddt = (0.5d0*dt2)/mass_colloid
	!write(*, fmt='(3F36.32)') dt2, ddt
	!call exit(0)
	do i=1,no_of_colloid
		dx=vel_colloid(3*i-2)*dt + f(3*i-2)*ddt
		dy=vel_colloid(3*i-1)*dt + f(3*i-1)*ddt
		dz=vel_colloid(3*i)*dt + f(3*i)*ddt

		pos_colloid(3*i-2)=pos_colloid(3*i-2) +dx
		pos_colloid(3*i-1)=pos_colloid(3*i-1) +dy
		pos_colloid(3*i)=pos_colloid(3*i) + dz


		if (pos_colloid(3*i-2).gt.llx) pos_colloid(3*i-2) = pos_colloid(3*i-2) - llx
		if (pos_colloid(3*i-1).gt.lly) pos_colloid(3*i-1) = pos_colloid(3*i-1) - lly
		if (pos_colloid(3*i).gt.llz) pos_colloid(3*i) = pos_colloid(3*i) - llz
		if (pos_colloid(3*i-2).lt.0.0d0 )   pos_colloid(3*i-2) = pos_colloid(3*i-2) +llx
        if (pos_colloid(3*i-1).lt.0.0d0 )   pos_colloid(3*i-1) = pos_colloid(3*i-1) +lly
       	if (pos_colloid(3*i).lt.0.0d0 )   pos_colloid(3*i) = pos_colloid(3*i) +llz
       	
    enddo
end subroutine update_pos_md
