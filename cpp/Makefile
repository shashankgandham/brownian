COMPILER = g++
EXECUTABLE = a.out
FLAGS = -O1 -pg
# Object Files for build
OBJS = .fluid_colloid_collision.o .initialize.o .neighbour_list.o .update.o .rotation_mpcd.o .run.o

$(EXECUTABLE) : active_brownian.cpp parameters.hpp $(OBJS)
	$(COMPILER) $(FLAGS) -o $(EXECUTABLE) parameters.hpp active_brownian.cpp $(OBJS)

# Object dependencies and compilation
.fluid_colloid_collision.o : fluid_colloid_collision.cpp
	 $(COMPILER) $(FLAGS) -c fluid_colloid_collision.cpp -o .fluid_colloid_collision.o

.initialize.o : initialize.cpp
	$(COMPILER) $(FLAGS) -c initialize.cpp -o .initialize.o

.neighbour_list.o : neighbour_list.cpp
	$(COMPILER) $(FLAGS) -c neighbour_list.cpp -o .neighbour_list.o

.update.o : update.cpp
	 $(COMPILER) $(FLAGS) -c update.cpp -o .update.o

.rotation_mpcd.o : rotation_mpcd.cpp
	 $(COMPILER) $(FLAGS) -c rotation_mpcd.cpp -o .rotation_mpcd.o

.run.o : run.cpp
	$(COMPILER) $(FLAGS) -c run.cpp -o .run.o

clean :
	rm -f .*.o 
