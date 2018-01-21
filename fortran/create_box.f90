subroutine create_box
	use all_parameters
	implicit none
	integer :: alx,aly,alz,i,j,k,l,m,mn,ii,pp,kk,temp,box
	do l=1,lz
		do m=1,ly
			do mn=1,lx
				box = (l-1)*lx*ly + (m-1)*lx + mn
				nbox = 0
				do kk = l,l+6
					k=kk-3
					if(k<=0) k= k + lz
					if(k>lz) k= k - lz
                	do pp=m,m+6
                		j=pp-3
                    	if(j<=0) j= j+ly
                    	if(j>ly) j=j-ly
                    	do ii=mn,mn+6
                    	   i=ii-3
                    	   if(i<=0) i= i + lx
                    	   if(i>lx) i= i - lx
                    	   temp = (k-1)*lx*ly + (j-1)*lx + i
                         
                         if(temp /= box ) then
                    	     nbox = nbox +1
                   			 box_neigh(nbox,box) = temp
                           endif
                   		enddo
               		enddo
          		enddo
    	   	enddo
		enddo
	enddo
end subroutine create_box
