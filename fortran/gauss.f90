subroutine gauss(v1,v2)   !! POLAR FORM OF BOX MULLER TRANSFORM
	use all_parameters
	implicit none
	real*8 :: x1,x2,z,v1,v2,ran1,sd
	z = dfloat(2)

	do while(z.gt.1.0d0)
		x1=2.0d0*ran1(zzzz) - 1
		x2=2.0d0*ran1(zzzz) - 1
		z=x1*x1+x2*x2
	end do
	z = dsqrt((-2.0d0*dlog(z))/z)

	v1 = x1*z*dsqrt(kbt/mass_fl) !gaussian with variance kbt/m and mean zero
	v2 = x2*z*dsqrt(kbt/mass_fl) !gaussian with variance kbt/m and mean zero
end subroutine gauss
