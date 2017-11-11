subroutine compute_force_md
	use all_parameters
	implicit none
	integer i,j,m
	real*8 x1,x2,y1,y2,z1,z2,x12,y12,z12,r,r1
	real*8 :: ffx,ffy,ffz, mag_f,r2,pot1

	potential_colloid =0.0d0
	f=0.0d0
	do i=1,no_of_colloid
		do j=1,n_neighbour(i)
			m = neighbour(j,i)
			!write(2,*) i,m
			x1=pos_colloid(3*i-2); x2=pos_colloid(3*m-2)
			y1=pos_colloid(3*i-1); y2=pos_colloid(3*m-1)
			z1=pos_colloid(3*i); z2=pos_colloid(3*m)

			x12=x1-x2; y12=y1-y2; z12=z1-z2
			x12 = x12 - llx*anint(x12*inv_llx)!minimum image convention
			y12 = y12 - lly*anint(y12*inv_lly)
			z12 = z12 - llz*anint(z12*inv_llz)
!           if(abs(x12).gt.llxby2) x12=(llx-abs(x12))*(x2-x1)/abs(x1-x2)
!           if(abs(y12).gt.llyby2) y12=(lly-abs(y12))*(y2-y1)/abs(y1-y2)
!           if(abs(z12).gt.llzby2) z12=(llz-abs(z12))*(z2-z1)/abs(z1-z2)
			r2 = x12*x12 + y12*y12 + z12*z12 ! NOTE
			if(r2.lt.r_cutoff2) then ! NOTE
				r = dsqrt(r2); r1=sig_colloid/r

				pot1 = (4.0d0*eps*(r1**12-r1**6) - ufc + fc*r)
				potential_colloid = potential_colloid + pot1
!				if(abs(pot1)>100.0d0) write(*,*)i,m,r
				mag_f = 4.0d0*eps*(12.0d0*sig_colloid12/r**13 - 6.0d0*sig_colloid6/r**7)  -fc !NOTE
				ffx = mag_f*x12/r; ffy = mag_f*y12/r; ffz = mag_f*z12/r

				f(3*i-2)=f(3*i-2) + ffx; f(3*i-1) = f(3*i-1) + ffy; f(3*i) = f(3*i) + ffz
				f(3*m-2)=f(3*m-2) - ffx; f(3*m-1) = f(3*m-1) - ffy; f(3*m) = f(3*m) - ffz
			endif
		enddo
	enddo
end subroutine compute_force_md

