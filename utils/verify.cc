#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <mpi.h>
#include <functional>
#include <ConfigFile.h>


extern "C" {
#include <libobjects.h>
#include <aoi_functions.h>
}

#define LINE_BUF_SIZE 256

#define EPS 0.00001

static std::string folder_name("./data");

static size_t clovis_block_size = 4096;

char clovis_local_addr[LINE_BUF_SIZE];
char clovis_ha_addr[LINE_BUF_SIZE];
char clovis_prof[LINE_BUF_SIZE];
char clovis_proc_fid[LINE_BUF_SIZE];
unsigned int tier = 0;

void load_config(std::string &rc_filename, int world_rank) {
	ConfigFile config(rc_filename, "=");
	std::string param_name;

	param_name = "LOCAL_ENDPOINT_ADDR" + std::to_string(world_rank);
	snprintf(clovis_local_addr, LINE_BUF_SIZE, "%s", config.read<std::string>(param_name).c_str());

	param_name = "LOCAL_PROC_FID" + std::to_string(world_rank);
	snprintf(clovis_proc_fid, LINE_BUF_SIZE, "%s", config.read<std::string>(param_name).c_str());

	snprintf(clovis_ha_addr, LINE_BUF_SIZE, "%s", config.read<std::string>("HA_ENDPOINT_ADDR").c_str());
	snprintf(clovis_prof, LINE_BUF_SIZE, "%s", config.read<std::string>("PROFILE_FID").c_str());

	printf("laddr: %s\t", clovis_local_addr);
	printf("ha_addr: %s\t", clovis_ha_addr);
	printf("prof: %s\t", clovis_prof);
	printf("proc fid: %s\n", clovis_proc_fid);
}

template<typename T>
int verify_array(T *truth, T *array, size_t size)
{
	for (size_t i = 0; i < size; i++) {
		if (std::abs(array[i] - truth[i]) > EPS) {
			std::cerr << i << " " << truth[i] << " vs " << array[i] << std::endl;
			return 1;
		}
	}
	return 0;
}

template<typename T>
void load_array(std::string &filename, T *array, size_t size)
{
	FILE *fp = fopen(filename.c_str(), "rb");
	assert (fp != NULL);
	memset(array, 0, sizeof(T) * size);
	fread(array, size, sizeof(T), fp);
	fclose(fp);
}

template<typename T>
void verify_particles(size_t species_count, size_t count, size_t iter, int world_rank, int world_size)
{
	T *data = nullptr;
	T *truth = (T*)malloc(sizeof(*truth) * count);
	std::string container_name;
	std::string object_name;
	std::string filename;
//	particle_species3-tracked_particles_5_v

	for (size_t species = 0; species < species_count; species++) {
		container_name = "particle_species" + std::to_string(species);
		Bucket *bkt = object_open_container(container_name.c_str());
		for (size_t cycle = 1; cycle <= iter; cycle++) {
			int get_rank, num_chunks, get_chunk_rank;
			int *get_dims, *get_chunk_dims; 
			size_t total_size;
	
			// x
			object_name = "tracked_particles_x_" + std::to_string(cycle);
			std::cout << "Preparing to get object: " << object_name << std::endl;
			data = (T*)object_get_chunk(bkt, object_name.c_str(), world_rank, FLOAT, &get_rank, &get_dims, &num_chunks, &get_chunk_rank, &get_chunk_dims);
	
			total_size = 1;
			for (size_t i = 0; i < get_chunk_rank; i++) total_size *= get_chunk_dims[i];

			std::cout << "Rank: " << get_rank << ": ";
			for (size_t i = 0; i < get_rank; i++) std::cout << get_dims[i] << " ";
			std::cout << std::endl;
	
			std::cout << "Num Chunks: " << num_chunks << " Chunk rank: " << get_chunk_rank << ": ";
			for (size_t i = 0; i < get_chunk_rank; i++) std::cout << get_chunk_dims[i] << " ";
			std::cout << std::endl;

			assert(get_rank == 1); assert(get_chunk_rank == 1);
			assert(get_dims[0] == count);
			assert(get_chunk_dims[0] == (count / world_size));
	
			filename = folder_name + "/particle_proc" + std::to_string(world_rank) + "_" + std::to_string(species) + "_x_" + std::to_string(cycle);
			std::cout << "load " << filename << std::endl;
			load_array<T>(filename, truth, count);
			assert(verify_array<T>(truth, data, total_size) == 0);
	
			free(get_dims); free(get_chunk_dims); free(data);
			std::cout << "Verified " << object_name.c_str() << " with " << filename.c_str() << std::endl;

			object_delete(bkt, object_name.c_str());
			std::cout << "Delete " << object_name << std::endl;

			// y
			object_name = "tracked_particles_y_" + std::to_string(cycle);
			std::cout << "Preparing to get object: " << object_name << std::endl;
			data = (T*)object_get_chunk(bkt, object_name.c_str(), world_rank, FLOAT, &get_rank, &get_dims, &num_chunks, &get_chunk_rank, &get_chunk_dims);
	
			total_size = 1;
			for (size_t i = 0; i < get_chunk_rank; i++) total_size *= get_chunk_dims[i];

			std::cout << "Rank: " << get_rank << ": ";
			for (size_t i = 0; i < get_rank; i++) std::cout << get_dims[i] << " ";
			std::cout << std::endl;
	
			std::cout << "Num Chunks: " << num_chunks << " Chunk rank: " << get_chunk_rank << ": ";
			for (size_t i = 0; i < get_chunk_rank; i++) std::cout << get_chunk_dims[i] << " ";
			std::cout << std::endl;

			assert(get_rank == 1); assert(get_chunk_rank == 1);
			assert(get_dims[0] == count);
			assert(get_chunk_dims[0] == (count / world_size));
	
			filename = folder_name + "/particle_proc" + std::to_string(world_rank) + "_" + std::to_string(species) + "_y_" + std::to_string(cycle);
			std::cout << "load " << filename << std::endl;
			load_array<T>(filename, truth, count);
			assert(verify_array<T>(truth, data, total_size) == 0);
	
			free(get_dims); free(get_chunk_dims); free(data);
			std::cout << "Verified " << object_name.c_str() << " with " << filename.c_str() << std::endl;

			object_delete(bkt, object_name.c_str());
			std::cout << "Delete " << object_name << std::endl;

			// z
			object_name = "tracked_particles_z_" + std::to_string(cycle);
			std::cout << "Preparing to get object: " << object_name << std::endl;
			data = (T*)object_get_chunk(bkt, object_name.c_str(), world_rank, FLOAT, &get_rank, &get_dims, &num_chunks, &get_chunk_rank, &get_chunk_dims);
	
			total_size = 1;
			for (size_t i = 0; i < get_chunk_rank; i++) total_size *= get_chunk_dims[i];

			std::cout << "Rank: " << get_rank << ": ";
			for (size_t i = 0; i < get_rank; i++) std::cout << get_dims[i] << " ";
			std::cout << std::endl;
	
			std::cout << "Num Chunks: " << num_chunks << " Chunk rank: " << get_chunk_rank << ": ";
			for (size_t i = 0; i < get_chunk_rank; i++) std::cout << get_chunk_dims[i] << " ";
			std::cout << std::endl;

			assert(get_rank == 1); assert(get_chunk_rank == 1);
			assert(get_dims[0] == count);
			assert(get_chunk_dims[0] == (count / world_size));
	
			filename = folder_name + "/particle_proc" + std::to_string(world_rank) + "_" + std::to_string(species) + "_z_" + std::to_string(cycle);
			std::cout << "load " << filename << std::endl;
			load_array<T>(filename, truth, count);
			assert(verify_array<T>(truth, data, total_size) == 0);
	
			free(get_dims); free(get_chunk_dims); free(data);
			std::cout << "Verified " << object_name.c_str() << " with " << filename.c_str() << std::endl;

			object_delete(bkt, object_name.c_str());
			std::cout << "Delete " << object_name << std::endl;

			// u
			object_name = "tracked_particles_u_" + std::to_string(cycle);
			std::cout << "Preparing to get object: " << object_name << std::endl;
			data = (T*)object_get_chunk(bkt, object_name.c_str(), world_rank, FLOAT, &get_rank, &get_dims, &num_chunks, &get_chunk_rank, &get_chunk_dims);
	
			total_size = 1;
			for (size_t i = 0; i < get_chunk_rank; i++) total_size *= get_chunk_dims[i];

			std::cout << "Rank: " << get_rank << ": ";
			for (size_t i = 0; i < get_rank; i++) std::cout << get_dims[i] << " ";
			std::cout << std::endl;
	
			std::cout << "Num Chunks: " << num_chunks << " Chunk rank: " << get_chunk_rank << ": ";
			for (size_t i = 0; i < get_chunk_rank; i++) std::cout << get_chunk_dims[i] << " ";
			std::cout << std::endl;

			assert(get_rank == 1); assert(get_chunk_rank == 1);
			assert(get_dims[0] == count);
			assert(get_chunk_dims[0] == (count / world_size));
	
			filename = folder_name + "/particle_proc" + std::to_string(world_rank) + "_" + std::to_string(species) + "_u_" + std::to_string(cycle);
			std::cout << "load " << filename << std::endl;
			load_array<T>(filename, truth, count);
			assert(verify_array<T>(truth, data, total_size) == 0);
	
			free(get_dims); free(get_chunk_dims); free(data);
			std::cout << "Verified " << object_name.c_str() << " with " << filename.c_str() << std::endl;

			object_delete(bkt, object_name.c_str());
			std::cout << "Delete " << object_name << std::endl;

			// v
			object_name = "tracked_particles_v_" + std::to_string(cycle);
			std::cout << "Preparing to get object: " << object_name << std::endl;
			data = (T*)object_get_chunk(bkt, object_name.c_str(), world_rank, FLOAT, &get_rank, &get_dims, &num_chunks, &get_chunk_rank, &get_chunk_dims);
	
			total_size = 1;
			for (size_t i = 0; i < get_chunk_rank; i++) total_size *= get_chunk_dims[i];

			std::cout << "Rank: " << get_rank << ": ";
			for (size_t i = 0; i < get_rank; i++) std::cout << get_dims[i] << " ";
			std::cout << std::endl;
	
			std::cout << "Num Chunks: " << num_chunks << " Chunk rank: " << get_chunk_rank << ": ";
			for (size_t i = 0; i < get_chunk_rank; i++) std::cout << get_chunk_dims[i] << " ";
			std::cout << std::endl;

			assert(get_rank == 1); assert(get_chunk_rank == 1);
			assert(get_dims[0] == count);
			assert(get_chunk_dims[0] == (count / world_size));
	
			filename = folder_name + "/particle_proc" + std::to_string(world_rank) + "_" + std::to_string(species) + "_v_" + std::to_string(cycle);
			std::cout << "load " << filename << std::endl;
			load_array<T>(filename, truth, count);
			assert(verify_array<T>(truth, data, total_size) == 0);
	
			free(get_dims); free(get_chunk_dims); free(data);
			std::cout << "Verified " << object_name.c_str() << " with " << filename.c_str() << std::endl;

			object_delete(bkt, object_name.c_str());
			std::cout << "Delete " << object_name << std::endl;

			// w
			object_name = "tracked_particles_w_" + std::to_string(cycle);
			std::cout << "Preparing to get object: " << object_name << std::endl;
			data = (T*)object_get_chunk(bkt, object_name.c_str(), world_rank, FLOAT, &get_rank, &get_dims, &num_chunks, &get_chunk_rank, &get_chunk_dims);
	
			total_size = 1;
			for (size_t i = 0; i < get_chunk_rank; i++) total_size *= get_chunk_dims[i];

			std::cout << "Rank: " << get_rank << ": ";
			for (size_t i = 0; i < get_rank; i++) std::cout << get_dims[i] << " ";
			std::cout << std::endl;
	
			std::cout << "Num Chunks: " << num_chunks << " Chunk rank: " << get_chunk_rank << ": ";
			for (size_t i = 0; i < get_chunk_rank; i++) std::cout << get_chunk_dims[i] << " ";
			std::cout << std::endl;

			assert(get_rank == 1); assert(get_chunk_rank == 1);
			assert(get_dims[0] == count);
			assert(get_chunk_dims[0] == (count / world_size));
	
			filename = folder_name + "/particle_proc" + std::to_string(world_rank) + "_" + std::to_string(species) + "_w_" + std::to_string(cycle);
			std::cout << "load " << filename << std::endl;
			load_array<T>(filename, truth, count);
			assert(verify_array<T>(truth, data, total_size) == 0);
	
			free(get_dims); free(get_chunk_dims); free(data);
			std::cout << "Verified " << object_name.c_str() << " with " << filename.c_str() << std::endl;

			object_delete(bkt, object_name.c_str());
			std::cout << "Delete " << object_name << std::endl;
		}
		object_close_container(bkt);
	}
	free(truth);
}

template <typename T>
void verify_quantities(size_t iter) 
{
	FILE *fp = nullptr;
	T *data = nullptr;
	T *truth = nullptr;
	std::string quantity_name, object_name;

	for (int quantity = 0; quantity < 3; quantity++) {
		for (int cycle = 1; cycle <= iter; cycle++) {
			size_t verify_count, verify_size;
	
			switch (quantity) {
				case 0: quantity_name = "rhoe"; break;
				case 1: quantity_name = "rhoi"; break;
				case 2: quantity_name = "rho_net"; break;
				default: throw "Invalid quantity ID!";
			}

	    		object_name = quantity_name + "_" + std::to_string(cycle);
			std::cout << "Preparing to get object: " << object_name << std::endl;
			int rc = aoi_get(object_name.c_str(), (char**)&data, &verify_count, &verify_size);
			assert(rc == 0);
			std::cout << "Retrieved object " << object_name << std::endl;

			truth = (T*)malloc(verify_count * verify_size);
	    		object_name = folder_name + "/" + quantity_name + "_" + std::to_string(cycle);
			load_array<T>(object_name, truth, verify_count);
			verify_array<T>(truth, data, verify_count);
			free(truth); free(data);
			std::cout << "Verified " << object_name << "!" << std::endl;

	    		object_name = quantity_name + "_" + std::to_string(cycle);
			aoi_delete(object_name.c_str());
			std::cout << "Delete " << object_name << std::endl;
		}
	}
}

int main(int argc, char *argv[]) {
	if (argc != 4) {
		printf("%s [rc_filename] [no. particles] [iter]\n", argv[0]);
		return 1;
	}
	int world_rank, world_size;
	size_t iter, no_particles;
	int rc;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	std::string rc_filename = std::string(argv[1]);
	no_particles = std::atoi(argv[2]);
	iter = std::atoi(argv[3]);

	load_config(rc_filename, world_rank);
	rc = aoi_init(clovis_local_addr, clovis_ha_addr, clovis_prof,
			clovis_proc_fid, clovis_block_size, tier);
	assert(rc == 0);
	std::cout << "Connected to mero..." << std::endl;
	verify_particles<float>(4, no_particles, iter, world_rank, world_size);
	if (world_rank == 0)
		verify_quantities<float>(iter);
	aoi_fini();
	MPI_Finalize();
	return 0;
}
