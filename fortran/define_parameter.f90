module all_parameters
	implicit none
    integer,parameter:: n=10,niter=21000
    double precision,parameter:: kbt=1.0d0,kbt1=1.0d0,ndt =0.1d0
	double precision,parameter:: dt=ndt/dfloat(n),pi=3.14159265d0

	integer,parameter:: thermostat=1,step_thermo=2
	integer,parameter:: lx=30,ly=30,lz=30
    double precision,parameter:: llx=dfloat(lx), lly=dfloat(ly), dt2=dt*dt
	double precision,parameter:: llz=dfloat(lz), llxby2=llx/2.0d0
	double precision,parameter:: llyby2=lly/2.0d0, llzby2=llz/2.0d0
    double precision,parameter:: inv_llx=1.0d0/llx,inv_lly=1.0d0/lly
	double precision,parameter:: inv_llz=1.0d0/llz
    double precision          :: ang_ke_colloid1,ang_ke_colloid
    integer                   :: zzz=-4280145,zzzz=77777,nn,qq
    
    integer,parameter         :: nbin=300
    double precision          :: dv=0.010d0         
    integer                   :: mb_vel(nbin)

    double precision,parameter:: mass_fl=1.0d0
    integer,parameter         :: no_of_fluid=lx*ly*lz*10,maxpart=100
    double precision          :: pos_fl(3*no_of_fluid),vel_fl(3*no_of_fluid),ke_fluid
    double precision,parameter:: vscale_fluid=dsqrt(12.0d0*kbt/mass_fl)
    double precision          ::scale_fac_mpcd
        
    integer,parameter           ::no_of_colloid=1
    double precision,parameter  ::mass_colloid=654.10d0,vscale_colloid=dsqrt(12.0d0*kbT1/mass_colloid)
    double precision,parameter  ::sig_colloid=5.0d0,r_cutoff=2.0**(1.0d0/6.0d0)*sig_colloid
    double precision,parameter  ::eps=1.0d0,space_limit=1.3d0*sig_colloid
    double precision,parameter  ::fc=4.0d0*eps*(12.0d0*(sig_colloid**12/r_cutoff**13)-6.0d0*(sig_colloid**6/r_cutoff**7))
    double precision,parameter  ::ufc=4.0d0*eps*((sig_colloid/r_cutoff)**12-(sig_colloid/r_cutoff)**6)+fc*r_cutoff
    double precision            ::scale_fac,ang_scale_fac
    integer						neighbour(200,no_of_colloid),n_neighbour(no_of_colloid)
    double precision            :: f_x(no_of_colloid),f_y(no_of_colloid),f_z(no_of_colloid),potential_colloid,ke_colloid
    double precision            ::f(3*no_of_colloid),r_cutoff2=r_cutoff**2
    double precision            ::sig_colloid12=sig_colloid**12,sig_colloid6=sig_colloid**6
    double precision,parameter  ::sig_colloidby2_square=sig_colloid**2*0.25d0
    double precision            ::old_force(3*no_of_colloid) !to copy the old force for updating colloid velocity
    double precision            ::pos_colloid(3*no_of_colloid),vel_colloid(3*no_of_colloid)
    double precision            ::ang_vel_colloid(3*no_of_colloid)
    double precision            ::rfx,rfy,rfz,rcx,rcy,rcz,mom_x,mom_y,mom_z
    double precision,parameter  ::sigma=0.80d0*sig_colloid,sigma_square=sigma**2 !actual diameter of colloid
    double precision            ::rsx,rsy,rsz,x_pos,y_pos,z_pos
    double precision            ::dx,dy,dz,ke_colloid1
    double precision,parameter  ::I_colloid=0.4*mass_colloid*sigma**2*0.25  !!MOMENT OF INERTIA
    double precision,parameter  ::ang_vscale_colloid=dsqrt(12.0d0*kbT1/I_colloid)
    double precision, parameter :: neigh_cutoff=3.0d0*sig_colloid,neigh_cutoff2=neigh_cutoff*neigh_cutoff !
    double precision,parameter  ::length_cutoff=sigma*0.50d0+2   !2 is more than 100*v*dt
    double precision,parameter  ::length_cutoff2=length_cutoff**2   !2 is more than 100*v*dt
    double precision            ::u1,u2,u3
    double precision            ::n_neighbour_fl(no_of_fluid)
    double precision,parameter  ::v0=0.04d0 ! 0 for no tumble
    integer           ::no_neigh(no_of_colloid),neigh_fl(10000,no_of_colloid)
    integer           ::box_neigh(500,lx*ly*lz),nbox
    double precision  ::ra(3*no_of_colloid),old_ra(3*no_of_colloid),rra(3*no_of_colloid)
    integer::size_cluster(no_of_colloid),n_cluster
    integer::identify(no_of_colloid)
    double precision ::dist(no_of_colloid,no_of_colloid)
    integer,parameter:: file_colloid=0
end module
