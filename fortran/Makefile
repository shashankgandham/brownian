COMPILER = gfortran
EXECUTABLE = active_brownian
FLAGS =
# Object Files for build
OBJS = \
define_parameter.o \
compute_force_md.o \
create_box.o \
fluid_colloid_collision.o \
initialize_activity.o \
initialize_colloid.o \
initialize_fluid.o \
neighbour_list_md.o \
neighbour_list_mpcd.o \
tumble.o \
update_activity_direction.o \
update_pos_md.o \
update_pos_mpcd.o \
update_velocity_colloid.o \
ran1.o \
rotation_mpcd.o \
run.o \
updown_velocity.o \
stochastic_reflection.o \
 gauss.o

$(EXECUTABLE) : active_brownian.f90 $(OBJS)
	$(COMPILER) $(FLAGS) -o $(EXECUTABLE) active_brownian.f90 $(OBJS)

# Object dependencies and compilation
define_parameter.o : define_parameter.f90
	 $(COMPILER) $(FLAGS) -c define_parameter.f90

compute_force_md.o : compute_force_md.f90
	$(COMPILER) $(FLAGS) -c compute_force_md.f90

create_box.o : create_box.f90
	$(COMPILER) $(FLAGS) -c create_box.f90

fluid_colloid_collision.o : fluid_colloid_collision.f90
	 $(COMPILER) $(FLAGS) -c fluid_colloid_collision.f90

initialize_activity.o : initialize_activity.f90
	$(COMPILER) $(FLAGS) -c initialize_activity.f90

initialize_colloid.o : initialize_colloid.f90
	$(COMPILER) $(FLAGS) -c initialize_colloid.f90

initialize_fluid.o : initialize_fluid.f90
	$(COMPILER) $(FLAGS) -c initialize_fluid.f90

neighbour_list_md.o : neighbour_list_md.f90
	$(COMPILER) $(FLAGS) -c neighbour_list_md.f90

neighbour_list_mpcd.o : neighbour_list_mpcd.f90
	 $(COMPILER) $(FLAGS) -c neighbour_list_mpcd.f90

tumble.o : tumble.f90
	$(COMPILER) $(FLAGS) -c  tumble.f90

update_activity_direction.o : update_activity_direction.f90
	 $(COMPILER) $(FLAGS) -c  update_activity_direction.f90

update_pos_md.o : update_pos_md.f90
	 $(COMPILER) $(FLAGS) -c  update_pos_md.f90

update_pos_mpcd.o : update_pos_mpcd.f90
	 $(COMPILER) $(FLAGS) -c update_pos_mpcd.f90

update_velocity_colloid.o : update_velocity_colloid.f90
	$(COMPILER) $(FLAGS) -c update_velocity_colloid.f90

ran1.o : ran1.f90
	 $(COMPILER) $(FLAGS) -c  ran1.f90

rotation_mpcd.o : rotation_mpcd.f90
	 $(COMPILER) $(FLAGS) -c rotation_mpcd.f90

run.o : run.f90
	$(COMPILER) $(FLAGS) -c run.f90

updown_velocity.o : updown_velocity.f90
	$(COMPILER) $(FLAGS) -c  updown_velocity.f90

stochastic_reflection.o : stochastic_reflection.f90
	$(COMPILER) $(FLAGS) -c  stochastic_reflection.f90

gauss.o : gauss.f90
	$(COMPILER) $(FLAGS) -c  gauss.f90
clean :
	rm -f *.mod *.o fort* *.dat
