       subroutine stochastic_reflection
       use all_parameters
       implicit none
       double precision::nx,ny,nz,unx,uny,unz,xx,m_beta,random_e
	   double precision::uuu,utx,uty,txx,tyy,tzz,ttt,utz,ran1,val,v1,v2
       m_beta=mass_fl/kbt
       !find the normal unit vector from colloid surface(nx,ny,nz)!!!!!!!!!!!!!!!!!!!1
       nx=(rsx)/dsqrt(rsx**2+rsy**2+rsz**2)
       ny=(rsy)/dsqrt(rsx**2+rsy**2+rsz**2)
       nz=(rsz)/dsqrt(rsx**2+rsy**2+rsz**2)

       !now generate magnitude of this unit vector from rayleigh distribution!!!!!!!!!??? Why??

        xx=ran1(zzzz)
        random_e=(1-xx)*(1-xx)
        val=dsqrt(-dlog(random_e)/m_beta)

       unx=val*nx     !un is normal velocity component after fluid colloid collision
       uny=val*ny
       unz=val*nz

       !now find a gaussian tangent vector to this normal!!!!!!!!!!!!!!!!!!!!!!!
       !first find the direction of tangent

       txx=ran1(zzzz)*llx
       tyy=ran1(zzzz)*lly
       tzz=ran1(zzzz)*llz

       txx=txx-rfx
       tyy=tyy-rfy
       tzz=tzz-rfz
                txx =txx - llx*anint(txx*inv_llx)
                tyy =tyy - lly*anint(tyy*inv_lly)
                tzz =tzz - llz*anint(tzz*inv_llz)

       utx=uny*tzz-tyy*unz
       uty=unz*txx-unx*tzz
       utz=unx*tyy-txx*uny
	   uuu=dsqrt(utx*utx+uty*uty+utz*utz)

       !unit vector along tangent
       utx=utx/uuu
       uty=uty/uuu
       utz=utz/uuu


       call gauss(v1,v2) !magnitude of tangent

       utx=v1*utx !ut is tangent velocity component after fluid colloid collision
       uty=v1*uty
       utz=v1*utz

      u1=unx+utx !x component of fluid velocity after fluid colloid collision
      u2=uny+uty
      u3=unz+utz

       end subroutine stochastic_reflection
