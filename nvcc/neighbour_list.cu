#include "parameters.cuh"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

inline __device__ __host__ point cmod(point a, point b) {
	if(a.x <=  0) a.x += b.x;  if(a.y <=  0) a.y += b.y;
	if(a.z <=  0) a.z += b.z; if(a.x > b.x) a.x -= b.x;
	if(a.y > b.y) a.y -= b.y; if(a.z > b.z) a.z -= b.z;
	return a;
}

inline __device__ point xyz(int cell, point len) {
	int x = len.x, y = len.y, px, py, pz;
	px = (cell%(x*y))%x;
	if(!px) px = 30; cell -= px;
	py = (cell%(x*y))/x;
	if(!py) py = 30; cell -= py;
	pz = cell/(x*y);
	return point(px, py, pz);
}

__global__ void d_create_box(int **box_neigh, point len) {
	int tbox, box, diff = len.y*len.z + len.x, nbox;
	point jiter, temp, iter;
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i <= len.prod()) {
		iter = xyz(i + diff, len);
		nbox = 0, box = (iter - point(0, 1, 1)).cell(len);
		jiter = iter - point(3, 3, 3);
		for(int j = 1; j <= 343; j++) {
			tbox = (cmod(jiter, len) - point(0, 1, 1)).cell(len);
			if(tbox != box) box_neigh[++nbox][box] = tbox;
			jiter.next(iter + point(3, 3, 3), point(1, 1, 1), iter - point(3, 3, 3));
		}
	}
}

void create_box() {
	int thr = 256, blk = (len.prod() + thr - 1)/thr;
	nbox = 342;
	d_create_box<<<blk, thr>>>(box_neigh, len);
}

__global__ void d_neighbour_list_md(int **neighbour, int *n_neighbour, point *pos_colloid, int no_of_colloid, double sig_colloid, point len) {
	double neigh_cutoff = 3.0*sig_colloid;
	point temp;
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	if(i < no_of_colloid) {
		for(int j = i + 1; j <= no_of_colloid; j++) {
			temp = img(pos_colloid[i] - pos_colloid[j], len);  
			if((temp*temp).sum() < power(neigh_cutoff,2)) 
				neighbour[++n_neighbour[i]][i] = j;
		}
	}
}
void neighbour_list_md() {
	int thr = 256, blk = (no_of_colloid + thr - 1)/thr;
	cudaMemset(n_neighbour, 0, sizeof(int)*(no_of_colloid + 2));  
	d_neighbour_list_md<<<blk, thr>>>(neighbour, n_neighbour, pos_colloid, no_of_colloid, sig_colloid, len);
}

__global__ void d_boxpart(int **box_part, int *fluid_no, int no_of_fluid, point *pos_fl, point len)  {
	int box_no;
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	if(i <= no_of_fluid) {
		box_no = 1 + pos_fl[i].cell(len);
		box_part[box_no][atomicAdd(&(fluid_no[box_no]), 1) + 1] = i;
	}
}
__global__ void neighfl_sync(int **neigh_fl, int *no_neigh, int no_of_colloid) {
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	if(i <= no_of_colloid)
		thrust::sort(thrust::seq, neigh_fl[i] + 1, neigh_fl[i] + no_neigh[i] + 1);
}

__global__ void d_neighbour_list_mpcd(int **box_part, int *fluid_no, int **box_neigh, int **neigh_fl, int *no_neigh,
	point *pos_colloid, point *pos_fl, int no_of_fluid, int no_of_colloid, int nbox, point len) {
	int mm, cbox;
	int j = blockIdx.x*blockDim.x + threadIdx.x + 1;
	if(j <= no_of_colloid) {
		no_neigh[j] = 0;
		cbox = 1 + pos_colloid[j].cell(len);
		for(int k = 1; k <= nbox; k++) {
			mm = box_neigh[k][cbox];
			for(int i = 1; i <= fluid_no[mm]; i++) {
				neigh_fl[j][atomicAdd(&no_neigh[j], 1) + 1] = box_part[mm][i];
			}
		}
	}
}


void neighbour_list_mpcd() {
	int thr = 512, blk = (no_of_fluid + thr - 1)/thr;
	cudaMemset(fluid_no, 0, sizeof(int)*(len.prod() + 2));
	d_boxpart<<<blk, thr>>>(box_part, fluid_no, no_of_fluid, pos_fl, len);
	blk = (len.prod() + thr - 1)/thr;
	blk = (no_of_colloid + thr - 1)/thr;
	d_neighbour_list_mpcd<<<blk, thr>>>(box_part, fluid_no, box_neigh, neigh_fl, no_neigh, pos_colloid, pos_fl, no_of_fluid, no_of_colloid, nbox, len);
	neighfl_sync<<<blk, thr>>>(neigh_fl, no_neigh, no_of_colloid);
}
