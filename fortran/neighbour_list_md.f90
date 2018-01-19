subroutine neighbour_list_md
    use all_parameters
    implicit none
    integer i,j
    double precision x1,x2,y1,y2,z1,z2,r,x12,y12,z12,r2

    n_neighbour =0
    neighbour=0

    do i=1,no_of_colloid-1
        do j=i+1,no_of_colloid
            x1=pos_colloid(3*i-2); y1=pos_colloid(3*i-1); 
            z1=pos_colloid(3*i);   x2=pos_colloid(3*j-2) 
            y2=pos_colloid(3*j-1); z2=pos_colloid(3*j)
            x12=x1-x2; y12=y1-y2; z12=z1-z2

            x12 = x12 - llx*anint(x12*inv_llx)
            y12 = y12 - lly*anint(y12*inv_lly)
            z12 = z12 - llz*anint(z12*inv_llz)
            
            r2=x12**2+y12**2+z12**2
            if(r2.lt.neigh_cutoff2) then
                n_neighbour(i)=n_neighbour(i)+1
                neighbour(n_neighbour(i),i)=j
            endif
        enddo
    enddo
end subroutine neighbour_list_md
