

#################### Program ###################################


EXEC=STO_proj_AOM_overlap        

CFLAGS = -Wall -O3

SRC_DIR = .

IFLAGS = -I/$(SRC_DIR)

OBJS +=  \
        $(SRC_DIR)/AOM_frag_cubes.o \
	$(SRC_DIR)/resolve_pvecs_per_fragment.o \
	$(SRC_DIR)/GTO_overlap.o \
	$(SRC_DIR)/AOM_overlap_calculation.o \
	$(SRC_DIR)/MakeLocalCoordsSystems.o \
	$(SRC_DIR)/STO_GTO_projection.o \
	$(SRC_DIR)/STO_proj_AOM_overlap.o \
	$(SRC_DIR)/anIres.o \
	$(SRC_DIR)/calculate_overlap_S_matrix.o \
	$(SRC_DIR)/create_cube_file.o \
	$(SRC_DIR)/equal.o \
	$(SRC_DIR)/initialize_STO.o \
	$(SRC_DIR)/locate_MO.o \
	$(SRC_DIR)/mulliken_functions.o \
	$(SRC_DIR)/read_CP2K_GTOs.o \
	$(SRC_DIR)/read_basis.o \
	$(SRC_DIR)/read_config_file.o \
	$(SRC_DIR)/read_xyz_convert_Ang_to_Bohr.o \
	$(SRC_DIR)/resolve_atomic_Z.o \
	$(SRC_DIR)/resolve_pvecs.o \
	$(SRC_DIR)/resolve_unique_species.o \
	$(SRC_DIR)/tokenize2.o \
        $(SRC_DIR)/simplex_obj_f.o \
        $(SRC_DIR)/run_simplex.o \
        $(SRC_DIR)/AOM_simplex_obj_f.o \
        $(SRC_DIR)/write_state.o \

STD_LIBS = -lm -llapack -lblas -lgfortran
#################### Target #####################################

all: $(EXEC)

$(EXEC): $(OBJS)
	gcc $(CFLAGS) $(IFLAGS) $(OBJS) -o $@ $(STD_LIBS)	

clean: 
	rm -f $(OBJS) $(EXEC)

