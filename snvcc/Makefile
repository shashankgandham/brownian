COMPILER = nvcc
EXECUTABLE = a.out
FLAGS = -O3 -G -Xcompiler -rdynamic
# Object Files for build
OBJS = .fluid_colloid_collision.o .initialize.o .neighbour_list.o .update.o .rotation_mpcd.o .run.o

$(EXECUTABLE) : active_brownian.cu $(OBJS)
	$(COMPILER) $(FLAGS) -o $(EXECUTABLE) active_brownian.cu $(OBJS)

# Object dependencies and compilation
.fluid_colloid_collision.o : fluid_colloid_collision.cu
	 $(COMPILER) $(FLAGS) -c fluid_colloid_collision.cu -o .fluid_colloid_collision.o

.initialize.o : initialize.cu
	$(COMPILER) $(FLAGS) -c initialize.cu -o .initialize.o

.neighbour_list.o : neighbour_list.cu
	$(COMPILER) $(FLAGS) -c neighbour_list.cu -o .neighbour_list.o

.update.o : update.cu
	 $(COMPILER) $(FLAGS) -c update.cu -o .update.o

.rotation_mpcd.o : rotation_mpcd.cu
	 $(COMPILER) $(FLAGS) -c rotation_mpcd.cu -o .rotation_mpcd.o

.run.o : run.cu
	$(COMPILER) $(FLAGS) -c run.cu -o .run.o

clean :
	rm -f .*.o *.txt
