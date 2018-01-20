subroutine update_pos_mpcd
    use all_parameters
    implicit none
    integer::  i
    do i=1,no_of_fluid
        pos_fl(3*i-2)=pos_fl(3*i-2)+vel_fl(3*i-2)*dt
        pos_fl(3*i-1)=pos_fl(3*i-1)+vel_fl(3*i-1)*dt
        pos_fl(3*i)  =pos_fl(3*i)+vel_fl(3*i)*dt
        pos_fl(3*i-2) = pos_fl(3*i-2)-llx*anint((pos_fl(3*i-2)-llxby2)*inv_llx)
        pos_fl(3*i-1) = pos_fl(3*i-1)-lly*anint((pos_fl(3*i-1)-llyby2)*inv_lly)
        pos_fl(3*i) = pos_fl(3*i)-llz*anint((pos_fl(3*i)-llzby2)*inv_llz)
        !if (nn == 209) write(*, fmt='(3F36.32)') pos_fl(3*i - 2), pos_fl(3*i - 1), pos_fl(3*i)
    enddo
end subroutine update_pos_mpcd
