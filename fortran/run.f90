subroutine run
	use all_parameters
	implicit none
	integer::i,ii,j,jj,mm,cnt(no_of_colloid),nbr(7000,no_of_colloid)
	integer::up_cnt(no_of_colloid),up_nbr(7000,no_of_colloid)
	double precision ::cos_b1,cos_b2,cos_b3,sin_b1,sin_b2,sin_b3
	double precision ::b1,b2,b3,bb1,bb2,bb3,ran1,v01,v02,v03
	double precision ::m11,m12,m13,m21,m22,m23,m31,m32,m33
	double precision :: vector_x,vector_y,vector_z,vector,dot
	double precision ::velx1,vely1,velz1,velx2,vely2,velz2,velx12,vely12,velz12
	double precision :: vv, theta,change,ke_before,ke_after,xx,vv1,vv2
	double precision ::up_velx2,up_vely2,up_velz2,aaa,delx,dely,delz
	double precision ::dumx,dumy,dumz, temp

	delx=0.0d0;dely=0.0d0;delz=0.0d0
	do i=1,no_of_colloid
		dumx=vel_colloid(3*i-2)
		dumy=vel_colloid(3*i-1)
		dumz=vel_colloid(3*i)
                
        delx=ra(3*i-2)*v0
		dely=ra(3*i-1)*v0
		delz=ra(3*i)*v0

                vel_colloid(3*i-2)=vel_colloid(3*i-2)+delx !-v0*old_ra(3*i-2)
		vel_colloid(3*i-1)=vel_colloid(3*i-1)+dely !-v0*old_ra(3*i-1)
		vel_colloid(3*i)=vel_colloid(3*i)+delz !-v0*old_ra(3*i)

		
		cnt(i)=0;up_cnt(i)=0
		do j=1,no_neigh(i)
			jj=neigh_fl(j,i)
			vector_x=pos_fl(3*jj-2)-pos_colloid(3*i-2)
			vector_y=pos_fl(3*jj-1)-pos_colloid(3*i-1)
			vector_z=pos_fl(3*jj)-pos_colloid(3*i)

			vector_x =vector_x - llx*anint(vector_x*inv_llx)
            vector_y =vector_y - lly*anint(vector_y*inv_lly)
            vector_z =vector_z - llz*anint(vector_z*inv_llz)

			vector=vector_x*vector_x+vector_y*vector_y+vector_z*vector_z
			vector_x = vector_x*vel_colloid(3*i-2)
			vector_y = vector_y*vel_colloid(3*i-1)
			vector_z = vector_z*vel_colloid(3*i)
			dot      = vector_x+vector_y+vector_z

        	if(vector.le.(sigma*0.50d0+0.5)**2.and.dot.le.0.0d0) then
         	    cnt(i)=cnt(i)+1
         		nbr(cnt(i),i)=jj
	        endif
        enddo

       do ii=1,cnt(i)
		   mm=nbr(ii,i)
		   xx=vel_fl(3*mm-2)
		   temp = mass_colloid/(mass_fl*dfloat(cnt(i)))
           dumx = delx*temp; dumy = dely*temp; dumz = delz*temp
           vel_fl(3*mm-2)=vel_fl(3*mm-2)-dumx
		   vel_fl(3*mm-1)=vel_fl(3*mm-1)-dumy
		   vel_fl(3*mm)=vel_fl(3*mm)-dumz
        enddo
	enddo
end subroutine run
