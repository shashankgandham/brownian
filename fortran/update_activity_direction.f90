       subroutine update_activity_direction
       use all_parameters
       implicit none
       integer::i
       real*8 ::b1,b2,b3,bb1,bb2,bb3,raa(3*no_of_colloid)
       REAL*8 ::m11,m12,m13,m21,m22,m23,m31,m32,m33
      do i=1,no_of_colloid
       b1=ang_vel_colloid(3*i-2)*dt
       b2=ang_vel_colloid(3*i-1)*dt
       b3=ang_vel_colloid(3*i)*dt
      
!FIRST MAKE A 3D ROTATION MATRIX      
       m11=cos(b2)*cos(b3)
       m12=-cos(b2)*sin(b3)
       m13=sin(b2)
       m21=sin(b1)*sin(b2)*cos(b3)+cos(b1)*sin(b3)
       m22=-sin(b1)*sin(b2)*sin(b3)+cos(b1)*cos(b3)
       m23=-sin(b1)*cos(b2)
       m31=-cos(b1)*sin(b2)*cos(b3)+sin(b1)*sin(b3)
       m32=cos(b1)*sin(b2)*sin(b3)+sin(b1)*cos(b3)
       m33=cos(b1)*cos(b2)

       bb1=m11*(m22*m33-m32*m23)
       bb2=-m12*(m21*m33-m31*m23)
       bb3=m13*(m21*m32-m31*m22)
!       write(*,*) bb1+bb2+bb3
!UPDATE THE DIRECTION OF ACTIVITY
       raa(3*i-2)=ra(3*i-2)
       raa(3*i-1)=ra(3*i-1)
       raa(3*i)=ra(3*i)

       ra(3*i-2)=m11*raa(3*i-2)+m12*raa(3*i-1)+m13*raa(3*i)
       ra(3*i-1)=m21*raa(3*i-2)+m22*raa(3*i-1)+m23*raa(3*i)
       ra(3*i)=m31*raa(3*i-2)+m32*raa(3*i-1)+m33*raa(3*i)
!      write(999,*) ra(3*i-2)**2 +ra(3*i-1)**2 +ra(3*i)**2
      enddo  

       end subroutine update_activity_direction
