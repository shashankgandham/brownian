subroutine tumble
    use all_parameters
    implicit none
    integer::i,ii,j,jj,mm,cnt(no_of_colloid),nbr(7000,no_of_colloid)
    integer::up_cnt(no_of_colloid),up_nbr(7000,no_of_colloid)
    double precision ::cos_b1,cos_b2,cos_b3,sin_b1,sin_b2,sin_b3
    double precision ::b1,b2,b3,bb1,bb2,bb3,ran1,v01,v02,v03
    double precision ::m11,m12,m13,m21,m22,m23,m31,m32,m33
    double precision :: vector_x,vector_y,vector_z,vector,dot
    double precision ::velx1,vely1,velz1,velx2,vely2,velz2,velx12,vely12
    double precision :: vv, theta,change,ke_before,ke_after,xx,vv1,vv2
    double precision ::up_velx2,up_vely2,up_velz2,aaa,delx,dely,delz
    double precision ::dumx,dumy,dumz,velz12

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
        ra(3*i-2) = ra(3*i-2) - llx*anint(ra(3*i-2)*inv_llx)
        ra(3*i-1) = ra(3*i-1) - llx*anint(ra(3*i-1)*inv_llx)
        ra(3*i) = ra(3*i) - llx*anint(ra(3*i)*inv_llx)

        aaa=dsqrt(ra(3*i-2)**2+ra(3*i-1)**2+ra(3*i)**2)
        
        ra(3*i-2)=ra(3*i-2)/aaa
        ra(3*i-1)=ra(3*i-1)/aaa
        ra(3*i)=ra(3*i)/aaa
    enddo
end subroutine tumble
