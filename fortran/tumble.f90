        subroutine tumble
        use all_parameters
        implicit none
        integer::i,ii,j,jj,mm,cnt(no_of_colloid),nbr(7000,no_of_colloid)
        integer::up_cnt(no_of_colloid),up_nbr(7000,no_of_colloid)
       double precision ::cos_b1,cos_b2,cos_b3,sin_b1,sin_b2,sin_b3
       double precision ::b1,b2,b3,bb1,bb2,bb3,ran1,v01,v02,v03
       REAL*8 ::m11,m12,m13,m21,m22,m23,m31,m32,m33
       double precision :: vector_x,vector_y,vector_z,vector,dot
       double precision ::velx1,vely1,velz1,velx2,vely2,velz2,velx12,vely12,velz12
       double precision :: vv, theta,change,ke_before,ke_after,xx,vv1,vv2
       double precision ::up_velx2,up_vely2,up_velz2,aaa,delx,dely,delz
       double precision ::dumx,dumy,dumz



!The bacterium has to tumble first to get the direction of new velocity
!the tumble angle is random in (0,2pi)
!      if(tumble.ne.0)then

       delx=0.0d0;dely=0.0d0;delz=0.0d0
      do i=1,no_of_colloid

          dumx=vel_colloid(3*i-2)
          dumy=vel_colloid(3*i-1)
          dumz=vel_colloid(3*i)

         b1=ran1(zzzz)*llx
         b2=ran1(zzzz)*lly
         b3=ran1(zzzz)*llz

              ra(3*i-2)= pos_colloid(3*i-2)-b1
              ra(3*i-1)= pos_colloid(3*i-1)-b2
              ra(3*i)= pos_colloid(3*i)-b3

                  ra(3*i-2) = ra(3*i-2) - llx*anint(ra(3*i-2)*inv_llx)!minimum image convention
                  ra(3*i-1) = ra(3*i-1) - llx*anint(ra(3*i-1)*inv_llx)!minimum image convention
                  ra(3*i) = ra(3*i) - llx*anint(ra(3*i)*inv_llx)!minimum image convention

              aaa=dsqrt(ra(3*i-2)**2+ra(3*i-1)**2+ra(3*i)**2)

              ra(3*i-2)=ra(3*i-2)/aaa !RA IS THE UNIT VECTOR TOWARDS THE ACTIVE VELOCITY
              ra(3*i-1)=ra(3*i-1)/aaa
              ra(3*i)=ra(3*i)/aaa


!UPDATE THE DIRECTION OF TUMBLE




!           vel_colloid(3*i-2)=vel_colloid(3*i-2)*ra(3*i-2) !-v0*old_ra(3*i-2)  
!           vel_colloid(3*i-1)=vel_colloid(3*i-1)*ra(3*i-1) !-v0*old_ra(3*i-1)  
!           vel_colloid(3*i)=vel_colloid(3*i)*ra(3*i) !-v0*old_ra(3*i)  

!          delx=vel_colloid(3*i-2)-dumx
!          dely=vel_colloid(3*i-1)-dumy
!          delz=vel_colloid(3*i)-dumz

! the added momenta has to be taken out from the neighbouring fluid
! to conserve momentum
!       cnt(i)=0;up_cnt(i)=0

!       write(*,*) ra(3*i-2),ra(3*i-2),ra(3*i)
!       do j=1,no_neigh(i)
!           jj=neigh_fl(j,i)



!          vector_x=pos_fl(3*jj-2)-pos_colloid(3*i-2)
!          vector_y=pos_fl(3*jj-1)-pos_colloid(3*i-1)
!          vector_z=pos_fl(3*jj)-pos_colloid(3*i)

!          vector_x =vector_x - llx*anint(vector_x*inv_llx)
!          vector_y =vector_y - lly*anint(vector_y*inv_lly)
!          vector_z =vector_z - llz*anint(vector_z*inv_llz)

!         vector=vector_x*vector_x+vector_y*vector_y+vector_z*vector_z

!          dot=vector_x*vel_colloid(3*i-2)+vector_y*vel_colloid(3*i-1)+vector_z*vel_colloid(3*i)


!        if(vector.le.(sigma*0.50d0+1)**2) then

!         cnt(i)=cnt(i)+1   !particles in the lower half and inside
!         nbr(cnt(i),i)=jj  ! spherical shell betn radii 2 to 3
!        endif



!        enddo

!       write(*,*)i,cnt(i)


!       do ii=1,cnt(i)

!          mm=nbr(ii,i)
!          xx=vel_fl(3*mm-2)

!         vel_fl(3*mm-2)=vel_fl(3*mm-2)-delx*mass_colloid/(mass_fl*dfloat(cnt(i)))
!         vel_fl(3*mm-1)=vel_fl(3*mm-1)-dely*mass_colloid/(mass_fl*dfloat(cnt(i)))
!         vel_fl(3*mm)=vel_fl(3*mm)-delz*mass_colloid/(mass_fl*dfloat(cnt(i)))

!       enddo

      enddo ! no_of_colloid

        end subroutine tumble

