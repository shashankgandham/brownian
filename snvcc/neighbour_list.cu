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
    for(int i = 1; i < no_of_colloid; i++) {
        for(int j = i + 1; j <= no_of_colloid; j++) {
            d_img(&temp, pos_colloid[i] - pos_colloid[j], len);  
            if((temp*temp).sum() < power(neigh_cutoff,2)) 
                neighbour[++n_neighbour[i]][i] = j;
        }
    }
}
__global__ void d_neighbour_list_mpcd(int **box_part, int *fluid_no, int **box_neigh, int **neigh_fl, int *no_neigh,
                    point *pos_colloid, point *pos_fl, int no_of_fluid, int no_of_colloid, int nbox, point len) {
    int box_no, mm, cbox;
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

void neighbour_list_mpcd() {
    memset(fluid_no, 0, sizeof(int)*(len.prod() + 2));
    point *d_pos_colloid, *d_pos_fl;
    int *d_fluid_no, **d_box_part, **d_neigh_fl, **d_box_neigh, *d_no_neigh;
    int *h_neigh_fl[10005], *h_box_part[128], *h_box_neigh[512];
	cudaMalloc(&d_pos_colloid, (no_of_colloid + 2)*sizeof(point)); 
    cudaMalloc(&d_pos_fl, (no_of_fluid + 2)*sizeof(point));
    cudaMalloc(&d_fluid_no, (no_of_fluid + 2)*sizeof(point));
    cudaMalloc(&d_no_neigh, (no_of_colloid + 2)*sizeof(point));
    cudaMalloc(&d_neigh_fl, (10005)*sizeof(int *));
    cudaMalloc(&d_box_part, (maxpart + 2)*sizeof(int *));
    cudaMalloc(&d_box_neigh, (512)*sizeof(int *));
	for(int i = 0; i <= 10000; i++) {
		cudaMalloc(&h_neigh_fl[i], (no_of_colloid + 2)*sizeof(int));
        cudaMemcpy(h_neigh_fl[i], neigh_fl[i], (no_of_colloid + 2)*sizeof(int), cudaMemcpyHostToDevice);
        if(i <= maxpart) {
            cudaMalloc(&h_box_part[i], (len.prod() + 2)*sizeof(int));
            cudaMemcpy(h_box_part[i], box_part[i], (len.prod() + 2)*sizeof(int), cudaMemcpyHostToDevice);
        }
        if(i <= 500) {
            cudaMalloc(&h_box_neigh[i], (len.prod() + 2)*sizeof(int));
            cudaMemcpy(h_box_neigh[i], box_neigh[i], (len.prod() + 2)*sizeof(int), cudaMemcpyHostToDevice);
        }
    }
    cudaMemcpy(d_pos_colloid, pos_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice ); 
    cudaMemcpy(d_pos_fl, pos_fl, (no_of_fluid + 2)*sizeof(point), cudaMemcpyHostToDevice);
    cudaMemcpy(d_fluid_no, fluid_no, (no_of_fluid + 2)*sizeof(point), cudaMemcpyHostToDevice);
    cudaMemcpy(d_no_neigh, no_neigh, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
    cudaMemcpy(d_neigh_fl, h_neigh_fl, (10005)*sizeof(int *), cudaMemcpyHostToDevice);
    cudaMemcpy(d_box_part, h_box_part, (maxpart + 2)*sizeof(int *), cudaMemcpyHostToDevice);
    cudaMemcpy(d_box_neigh, h_box_neigh, (512)*sizeof(int *), cudaMemcpyHostToDevice);
    d_neighbour_list_mpcd<<<1, 1>>>(d_box_part, d_fluid_no, d_box_neigh, d_neigh_fl, d_no_neigh, d_pos_colloid, d_pos_fl, 
                            no_of_fluid, no_of_colloid, nbox, len);    
   for(int i = 0; i <= 10000; i++) {
        cudaMemcpy(neigh_fl[i], h_neigh_fl[i], (no_of_colloid + 2)*sizeof(int), cudaMemcpyDeviceToHost);
        cudaFree(neigh_fl[i]);
        if(i <= 500) cudaFree(h_box_neigh[i]);
        if(i <= maxpart) {    
            cudaMemcpy(box_part[i], h_box_part[i], (len.prod() + 2)*sizeof(int), cudaMemcpyDeviceToHost);
            cudaFree(h_box_part[i]);
        }
    }
    cudaMemcpy(fluid_no, d_fluid_no, (no_of_fluid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
    cudaMemcpy(no_neigh, d_no_neigh, (no_of_colloid + 2)*sizeof(point), cudaMemcpyDeviceToHost);
    cudaFree(d_fluid_no), cudaFree(d_no_neigh), cudaFree(d_pos_colloid), cudaFree(d_pos_fl);
    cudaFree(d_box_part), cudaFree(d_box_neigh), cudaFree(d_neigh_fl);
}

void neighbour_list_md() {
    memset(n_neighbour, 0, sizeof(int)*(no_of_colloid + 2));
    point *d_pos;
	int *d_n_neighbour, **d_neighbour, *h_neighbour[256];
	cudaMalloc(&d_n_neighbour, (no_of_colloid + 2)*sizeof(int)); 
    cudaMalloc(&d_neighbour, 256*sizeof(int *));
    cudaMalloc(&d_pos, sizeof(point)*(no_of_colloid + 2));
	for(int i = 0; i <= 200; i++) {
		cudaMalloc(&h_neighbour[i], (no_of_colloid + 2)*sizeof(int));
		cudaMemcpy(h_neighbour[i], neighbour[i], (no_of_colloid + 2)*sizeof(int), cudaMemcpyHostToDevice);
    }
    cudaMemcpy(d_neighbour, h_neighbour, 256*sizeof(int *), cudaMemcpyHostToDevice);
    cudaMemcpy(d_pos, pos_colloid, (no_of_colloid + 2)*sizeof(point), cudaMemcpyHostToDevice);
    cudaMemcpy(d_n_neighbour, n_neighbour, (no_of_colloid + 2)*sizeof(int), cudaMemcpyHostToDevice);  
  	d_neighbour_list_md<<<1, 1>>>(d_neighbour, d_n_neighbour, d_pos, no_of_colloid, sig_colloid, len);
    for(int i = 0; i <= 200; i++) {
		cudaMemcpy(neighbour[i], h_neighbour[i], (no_of_colloid + 2)*sizeof(int), cudaMemcpyDeviceToHost);
        cudaFree(h_neighbour[i]);
    }
    cudaFree(d_neighbour), cudaFree(d_n_neighbour), cudaFree(d_pos);
}