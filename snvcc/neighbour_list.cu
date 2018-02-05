#include "parameters.cuh"

inline point cmod(point a, point b) {
    if(a.x <=  0) a.x += b.x;  if(a.y <=  0) a.y += b.y;
    if(a.z <=  0) a.z += b.z; if(a.x > b.x) a.x -= b.x;
    if(a.y > b.y) a.y -= b.y; if(a.z > b.z) a.z -= b.z;
    return a;
}
void create_box() {
    int tbox, box;
    point jiter, iter = point(1, 1, 1), temp;
    for(int i = 1; i <= len.prod(); i++, iter.next(len)) {
        nbox = 0, box = (iter - point(0, 1, 1)).cell(len);
        jiter = iter - point(3, 3, 3);
        for(int j = 1; j <= 343; j++) {
            tbox = (cmod(jiter, len) - point(0, 1, 1)).cell(len);
            if(tbox != box) box_neigh[++nbox][box] = tbox;
            jiter.next(iter + point(3, 3, 3), point(1, 1, 1), iter - point(3, 3, 3));
        }
    }
}
__global__ void d_neighbour_list_md(int **neighbour, int *n_neighbour, point *pos_colloid, int no_of_colloid, double sig_colloid, point len) {
    double neigh_cutoff = 3.0*sig_colloid;
    point temp;
    memset(n_neighbour, 0, sizeof(int)*(no_of_colloid + 2));
    int i = blockIdx.x + 1, j = blockIdx.x + blockIdx.y + 1;
    if(j <= no_of_colloid) {
      d_img(&temp, pos_colloid[i] - pos_colloid[j], len);  
      if((temp*temp).sum() < pow(neigh_cutoff,2)) 
        neighbour[++n_neighbour[i]][i] = j;
    }
}
void neighbour_list_md() {
    point *d_pos;
    memset(n_neighbour, 0, sizeof(int)*(no_of_colloid + 2));
	int *d_n_neighbour, **d_neighbour, *h_neighbour[256];
	cudaMalloc(&d_n_neighbour, (no_of_colloid + 2)*sizeof(int)); 
    cudaMalloc(&d_neighbour, 256*sizeof(int *));
    cudaMalloc(&d_pos, sizeof(point)*(no_of_colloid + 2));
	for(int i = 0; i <= 200; i++) {
		cudaMalloc(&h_neighbour[i], (no_of_colloid + 2)*sizeof(int));
		cudaMemcpy(h_neighbour[i], neighbour[i], (no_of_colloid + 2)*sizeof(int), cudaMemcpyHostToDevice);
	}
  	cudaMemcpy(d_pos, pos_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
  	d_neighbour_list_md<<<no_of_colloid + 1, 1>>>(d_neighbour, d_n_neighbour, d_pos, no_of_colloid, sig_colloid, len);
  	cudaMemcpy(neighbour, d_neighbour, 256*(no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
}

void neighbour_list_mpcd() {
    int box_no, mm, cbox;
    memset(fluid_no, 0, sizeof(int)*(len.prod() + 2));
    for(int	i = 1; i <= no_of_fluid; i++) {
        box_no = 1 + pos_fl[i].cell(len);
        box_part[++fluid_no[box_no]][box_no] = i;
    }
    for(int j = 1; j <= no_of_colloid; j++) {
        no_neigh[j] = 0;
        cbox = 1 + pos_colloid[j].cell(len);
        for(int k = 1; k <= nbox; k++) {
            mm = box_neigh[k][cbox];
            for(int i = 1; i <= fluid_no[mm]; i++) {
                neigh_fl[++no_neigh[j]][j] = box_part[i][mm];
            }
        }
    }
}