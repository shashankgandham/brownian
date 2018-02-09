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
            temp = img(pos_colloid[i] - pos_colloid[j], len);  
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

void neighbour_list_md() {
    cudaMemset(n_neighbour, 0, sizeof(int)*(no_of_colloid + 2));  
    d_neighbour_list_md<<<1, 1>>>(neighbour, n_neighbour, pos_colloid, no_of_colloid, sig_colloid, len);
}

void neighbour_list_mpcd() {
    cudaMemset(fluid_no, 0, sizeof(int)*(len.prod() + 2));
    d_neighbour_list_mpcd<<<1, 1>>>(box_part, fluid_no, box_neigh, neigh_fl, no_neigh, pos_colloid, pos_fl, no_of_fluid, no_of_colloid, nbox, len);
}
