       subroutine updown_velocity
       use all_parameters
       implicit none
       double precision ::up_velx2,up_vely2,up_velz2
       double precision :: vector_x,vector_y,vector_z,vector,dot
       double precision ::velx2,vely2,velz2,velx12,vely12,velz12
        integer::i,ii,j,jj,mm,cnt(no_of_colloid),nbr(7000,no_of_colloid)
        integer::up_cnt(no_of_colloid),up_nbr(7000,no_of_colloid)

           do i=1,no_of_colloid

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

          dot=vector_x*vel_colloid(3*i-2)+vector_y*vel_colloid(3*i-1)+vector_z*vel_colloid(3*i)
       !   dot=dot/(dsqrt(vector)*dsqrt(vel_colloid(3*i-2)**2+vel_colloid(3*i-1)**2+vel_colloid(3*i)**2))
       !   dot=(dot*180.0d0)/acos(-1.0d0)

        if(vector.le.(sigma*0.50d0+0.5)**2.and.dot.le.0.0d0) then

         cnt(i)=cnt(i)+1   !particles in the lower half and inside
         nbr(cnt(i),i)=jj  ! spherical shell betn radii 2 to 3
        endif

        if(vector.le.(sigma*0.50d0+0.1)**2.and.dot.ge.0.0d0) then

        up_cnt(i)=up_cnt(i)+1   !particles in the upper half and inside
        up_nbr(up_cnt(i),i)=jj  ! spherical shell betn radii 2 to 3

        endif

        enddo
       velx2=0.0d0
       vely2=0.0d0
       velz2=0.0d0

      if (cnt(i)>0.and.up_cnt(i)>0) then
       do ii=1,cnt(i)

          mm=nbr(ii,i)


         velx2=velx2+vel_fl(3*mm-2)
         vely2=vely2+vel_fl(3*mm-1)
         velz2=velz2+vel_fl(3*mm)
       enddo

       up_velx2=0.0d0
       up_vely2=0.0d0
       up_velz2=0.0d0

       do ii=1,up_cnt(i)
          mm=up_nbr(ii,i)
         up_velx2=up_velx2+vel_fl(3*mm-2)
         up_vely2=up_vely2+vel_fl(3*mm-1)
         up_velz2=up_velz2+vel_fl(3*mm)
       enddo

       up_velx2=up_velx2/dfloat(up_cnt(i)) !avg vel along colloid velocity
       up_vely2=up_vely2/dfloat(up_cnt(i)) 
       up_velz2=up_velz2/dfloat(up_cnt(i)) 

        velx2=velx2/dfloat(cnt(i))     !avg vel opposite to colloid velocity
        vely2=vely2/dfloat(cnt(i))
        velz2=velz2/dfloat(cnt(i))

       
       up_velx2=up_velx2-vel_colloid(3*i-2) !relative avg vel along colloid velocity
       up_vely2=up_vely2-vel_colloid(3*i-1) 
       up_velz2=up_velz2-vel_colloid(3*i) 

        velx2=velx2-vel_colloid(3*i-2)     !relative avg vel opposite to colloid velocity
        vely2=vely2-vel_colloid(3*i-1)
        velz2=velz2-vel_colloid(3*i)

      endif
      enddo
       end subroutine 
