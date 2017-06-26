      subroutine update_pos_mpcd    !position update for all fluid particles 
                                    !for md time step
      use all_parameters
      implicit none
      integer::  i

       do i=1,no_of_fluid
            pos_fl(3*i-2)=pos_fl(3*i-2)+vel_fl(3*i-2)*dt
            pos_fl(3*i-1)=pos_fl(3*i-1)+vel_fl(3*i-1)*dt
            pos_fl(3*i)  =pos_fl(3*i)+vel_fl(3*i)*dt

         pos_fl(3*i-2) = pos_fl(3*i-2) - llx*anint((pos_fl(3*i-2)-llxby2)*inv_llx)!pbc using minimum image
         pos_fl(3*i-1) = pos_fl(3*i-1) - lly*anint((pos_fl(3*i-1)-llyby2)*inv_lly)
         pos_fl(3*i) = pos_fl(3*i) - llz*anint((pos_fl(3*i)-llzby2)*inv_llz)


!          if(pos_fl(3*i-2).gt.llx)   pos_fl(3*i-2)=pos_fl(3*i-2)-llx    ! PBC
!          if(pos_fl(3*i-2).le.0.0d0) pos_fl(3*i-2)=pos_fl(3*i-2)+llx
!          if(pos_fl(3*i-1).gt.lly)   pos_fl(3*i-1)=pos_fl(3*i-1)-lly
!          if(pos_fl(3*i-1).le.0.0d0) pos_fl(3*i-1)=pos_fl(3*i-1)+lly
!          if(pos_fl(3*i).gt.llz)     pos_fl(3*i)=pos_fl(3*i)-llz
!          if(pos_fl(3*i).le.0.0d0)   pos_fl(3*i)=pos_fl(3*i)+llz

       enddo

      end subroutine update_pos_mpcd
