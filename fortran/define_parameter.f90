! generate data to find msd

       module all_parameters
        implicit none
!SIMULATION  SPECIFICATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer,parameter ::n=10,niter=21000
        real*8,parameter  ::kbt=1.0d0,kbt1=1.0d0,ndt =0.1d0,dt=ndt/dfloat(n),pi=3.14159265d0
! kbt is system temp after equilibrium.kbt1 is to initiate the process at different temp

        integer,parameter :: thermostat=1,step_thermo=2  !stepthermo in terms of mpc step
        integer,parameter ::lx=30,ly=30,lz=30,dt2=dt*dt
        real*8,parameter  :: llx=dfloat(lx), lly=dfloat(ly), llz=dfloat(lz)
        real*8,parameter  :: llxby2=llx/2.0d0, llyby2=lly/2.0d0, llzby2=llz/2.0d0
        real*8,parameter  :: inv_llx=1.0d0/llx,inv_lly=1.0d0/lly,inv_llz=1.0d0/llz
        real*8            ::ang_ke_colloid1,ang_ke_colloid
        integer           :: zzz=-4280145,zzzz=77777,nn,qq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer,parameter ::nbin=300 !for maxwell velocity
        real*8            ::dv=0.010d0                     !distribution
        integer           ::mb_vel(nbin)
!FLUID SPECIFICATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        real*8,parameter ::mass_fl=1.0d0
        integer,parameter::no_of_fluid=lx*ly*lz*10,maxpart=100 !maximum no of fluid particle in one mpc cell
        real*8           ::pos_fl(3*no_of_fluid),vel_fl(3*no_of_fluid),ke_fluid
        real*8,parameter ::vscale_fluid=dsqrt(12.0d0*kbt/mass_fl)
        real*8           ::scale_fac_mpcd
!COLLOID SPECIFICATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        integer,parameter ::no_of_colloid=1
        real*8,parameter  ::mass_colloid=654.10d0,vscale_colloid=dsqrt(12.0d0*kbT1/mass_colloid)
        real*8,parameter  ::sig_colloid=5.0d0,r_cutoff=2.0**(1.0d0/6.0d0)*sig_colloid,eps=1.0d0,space_limit=1.3d0*sig_colloid
        real*8,parameter  ::fc=4.0d0*eps*(12.0d0*(sig_colloid**12/r_cutoff**13)-6.0d0*(sig_colloid**6/r_cutoff**7))
        real*8,parameter  ::ufc=4.0d0*eps*((sig_colloid/r_cutoff)**12-(sig_colloid/r_cutoff)**6)+fc*r_cutoff
        real*8            ::scale_fac,ang_scale_fac
        integer           :: neighbour(200,no_of_colloid),n_neighbour(no_of_colloid)
        real*8            :: f_x(no_of_colloid),f_y(no_of_colloid),f_z(no_of_colloid),potential_colloid,ke_colloid
        real*8            ::f(3*no_of_colloid),r_cutoff2=r_cutoff**2
        real*8            ::sig_colloid12=sig_colloid**12,sig_colloid6=sig_colloid**6
        real*8,parameter  ::sig_colloidby2_square=sig_colloid**2*0.25d0
        real*8            ::old_force(3*no_of_colloid) !to copy the old force for updating colloid velocity
        real*8            ::pos_colloid(3*no_of_colloid),vel_colloid(3*no_of_colloid)
        real*8            ::ang_vel_colloid(3*no_of_colloid)
        real*8            ::rfx,rfy,rfz,rcx,rcy,rcz,mom_x,mom_y,mom_z
        real*8,parameter  ::sigma=0.80d0*sig_colloid,sigma_square=sigma**2 !actual diameter of colloid
        real*8            ::rsx,rsy,rsz,x_pos,y_pos,z_pos
        real*8            ::dx,dy,dz,ke_colloid1
!sig_colloid is minimum distance between two colloids,r_cutoff is range of L-J force, space_limit is initial separation
!between any two colloids!  NOTE: SPACE_LIMIT MUST BE LESS THAN R_CUTOFF
!eps is epsilon in L-J.

        real*8,parameter  ::I_colloid=0.4*mass_colloid*sigma**2*0.25  !!MOMENT OF INERTIA
        real*8,parameter  ::ang_vscale_colloid=dsqrt(12.0d0*kbT1/I_colloid)
        real*8, parameter :: neigh_cutoff=3.0d0*sig_colloid,neigh_cutoff2=neigh_cutoff*neigh_cutoff !
!neigh_cutoff is the cutoff for constructing neighbour list
        real*8,parameter  ::length_cutoff=sigma*0.50d0+2   !2 is more than 100*v*dt
        real*8,parameter  ::length_cutoff2=length_cutoff**2   !2 is more than 100*v*dt
        real*8            ::u1,u2,u3
        real*8            ::n_neighbour_fl(no_of_fluid)
        integer           ::no_neigh(no_of_colloid),neigh_fl(10000,no_of_colloid)
        integer           ::box_neigh(500,lx*ly*lz),nbox
!FOR ACTIVITY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8            ::ra(3*no_of_colloid),old_ra(3*no_of_colloid),rra(3*no_of_colloid)
        real*8,parameter  ::v0=0.04d0 ! 0 for no tumble
!FOR CLUSTER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer::size_cluster(no_of_colloid),n_cluster
        integer::identify(no_of_colloid)
        real*8 ::dist(no_of_colloid,no_of_colloid)
        integer,parameter:: file_colloid=0 !set 1 to read initial
                                           !colloid and fluid from file

       end module

