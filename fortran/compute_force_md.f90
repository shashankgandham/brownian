subroutine compute_force_md
	use all_parameters
	implicit none
	integer i,j,m,l
	double precision x1,x2,y1,y2,z1,z2,x12,y12,z12,r,r1
    double precision :: ffx,ffy,ffz, mag_f,pot1, r2
    double precision t1, t2, t3, r13, r7
    potential_colloid =0.0d0
    f=0.0d0
	do i=1,no_of_colloid
		do j=1,n_neighbour(i)
			m = neighbour(j,i)
			x1=pos_colloid(3*i-2); x2=pos_colloid(3*m-2)
            y1=pos_colloid(3*i-1); y2=pos_colloid(3*m-1)
            z1=pos_colloid(3*i); z2=pos_colloid(3*m)

            x12=x1-x2; y12=y1-y2; z12=z1-z2
            x12 = x12 - llx*anint(x12*inv_llx)!minimum image convention
            y12 = y12 - lly*anint(y12*inv_lly)
            z12 = z12 - llz*anint(z12*inv_llz)
            r2=x12*x12+y12*y12+z12*z12 ! NOTE

            if(r2.lt.r_cutoff2) then ! NOTE-
	            r=dsqrt(r2); r1=sig_colloid/r
                r13 = 1.0d0; r7 = 1.0d0
                do l = 1, 13
                    r13 = r13*r
                enddo
                do l = 1, 7
                    r7 = r7*r
                enddo
                pot1=(4.0d0*eps*(r1**12-r1**6) -ufc + fc*r)
                potential_colloid=potential_colloid+pot1
                t1 = sig_colloid12/r13; t2 = sig_colloid6/r7
                mag_f = 4.0d0*eps*(12.0d0*t1 - 6.0d0*t2) - fc !NOTE
                ffx = (x12*mag_f)/r; ffy = (y12*mag_f)/r; ffz = (z12*mag_f)/r
                f(3*i-2)=f(3*i-2) + ffx; f(3*i-1) = f(3*i-1) + ffy; f(3*i) = f(3*i) + ffz
                f(3*m-2)=f(3*m-2) - ffx; f(3*m-1) = f(3*m-1) - ffy; f(3*m) = f(3*m) - ffz
            endif
		enddo
	enddo
end subroutine compute_force_md
