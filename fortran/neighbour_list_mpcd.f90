subroutine neighbour_list_mpcd
    use all_parameters
    implicit none
    integer::i,j,k,ii,box_no,fluid_no(lx*ly*lz)
    integer::p,l,lxly,mm, box_part(maxpart, lx*ly*lz)
    integer::bx1,bx2,by1,by2,bz1,bz2
    double precision ::x1,x2,y1,y2,z1,z2,sigmaby2
    integer::di,idx,idy,idz,int_sigma
    integer :: x1c,y1c,z1c,cbox

    sigmaby2=sigma*0.50d0
    fluid_no=0; box_part=0
    lxly= lx*ly
    do i=1,no_of_fluid
        box_no = 1+int(pos_fl(3*i-2))+lx*int(pos_fl(3*i-1))+lxly*int(pos_fl(3*i))
        fluid_no(box_no) = fluid_no(box_no) + 1
        j = fluid_no(box_no)
        box_part(j,box_no) = i

    enddo
    do j=1,no_of_colloid
        no_neigh(j)=0
        x1=int(pos_colloid(3*j-2))
        y1=int(pos_colloid(3*j-1))
        z1=int(pos_colloid(3*j))
        cbox=1+x1+y1*lx+z1*lx*ly
        do k=1,nbox
            mm = box_neigh(k,cbox)
            do i=1,fluid_no(mm)
                ii=box_part(i,mm)
                no_neigh(j)=no_neigh(j)+1
                neigh_fl(no_neigh(j),j)=ii
                
            enddo
        enddo
    enddo
end subroutine neighbour_list_mpcd
