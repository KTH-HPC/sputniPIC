#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <mpi.h>
#include <functional>


extern "C" {
#include <libobjects.h>
#include <aoi_functions.h>
}

#define LINE_BUF_SIZE 256

#define EPS 0.00001

static size_t clovis_block_size = 4096;

char clovis_local_addr[LINE_BUF_SIZE];
char clovis_ha_addr[LINE_BUF_SIZE];
char clovis_prof[LINE_BUF_SIZE];
char clovis_proc_fid[LINE_BUF_SIZE];
char rc_filename[] = "sagerc";
unsigned int tier = 0;

int get_line(FILE *fp, char *dst) {
	char buffer[LINE_BUF_SIZE] = {0};
	while (1) {
		assert(fgets(buffer, LINE_BUF_SIZE, fp) != NULL);
		if (strlen(buffer) == 0 || buffer[0] == '#' || buffer[0] == '\n')
			continue;
		char *end = buffer + strlen(buffer) - 1;
		while (*end == '\n')
			*(end--) = '\0';
		strcpy(dst, buffer);
		break;
	}
	return 1;
}

void load_config() {
	FILE *rc_file = fopen(rc_filename, "r");
	assert(rc_file != NULL);
	get_line(rc_file, clovis_local_addr);
	printf("laddr: %s\t", clovis_local_addr);
	get_line(rc_file, clovis_ha_addr);
	printf("ha_addr: %s\t", clovis_ha_addr);
	get_line(rc_file, clovis_prof);
	printf("prof: %s\t", clovis_prof);
	get_line(rc_file, clovis_proc_fid);
	printf("proc fid: %s\n", clovis_proc_fid);
	fclose(rc_file);
}

void verify_particles(size_t species_count, size_t count, size_t iter)
{
	FILE *fp = nullptr;
	float *data = nullptr;
	data = (float*)malloc(count * sizeof(float));
	size_t verify_count, verify_size;
	float *aoi_verify_data;
//	tracked_particles_proc0_12_v2
	for (size_t species = 0; species < species_count; species++) {
		for (size_t cycle = 1; cycle <= iter; cycle++) {
			std::string object_name = "data/tracked_particles_proc" + std::to_string(0) + "_" +
        		          std::to_string(cycle) + "_" + "x" + std::to_string(species);
			std::cout << "Try to verify " << object_name << "..." << std::endl;
			fp = fopen(object_name.c_str(), "rb");
			memset(data, 0, sizeof(float) * count);
			fread(data, count, sizeof(float), fp);
			fclose(fp);
			int rc = aoi_get(object_name.c_str(), (char**)&aoi_verify_data, &verify_count, &verify_size);
			assert(rc == 0);
			assert(verify_count == count);
			assert(verify_size == sizeof(float));
			for (size_t i = 0; i < count; i++) {
				// fabs
				if (std::abs(data[i] - aoi_verify_data[i]) > EPS) {
					std::cerr << i << " " << data[i] << " " << aoi_verify_data[i] << std::endl;
				}
			}
			free((char*)aoi_verify_data);
			std::cout << "Verified " << object_name.c_str() << std::endl;
		}
	}
	free(data);
}

void verify_quantities(size_t iter) 
{
	FILE *fp = nullptr;
	float *data = nullptr;
	float *buffer = nullptr;
	bucket *bkt = object_open_container("quantities");
	std::string quantity_name, filename;

	for (int quantity = 0; quantity < 3; quantity++) {
		for (int cycle = 1; cycle <= iter; cycle++) {
			int get_rank, num_chunks, get_chunk_rank;
			int *get_dims, *get_chunk_dims; 
	
			switch (quantity) {
				case 0: quantity_name = "rhoe"; break;
				case 1: quantity_name = "rhoi"; break;
				case 2: quantity_name = "rho_net"; break;
				default: throw "Invalid quantity ID!";
			}

	    		filename = quantity_name + "_" + std::to_string(cycle);
			std::cout << "Preparing to get object: " << filename << std::endl;
			data = (float*)object_get_chunk(bkt, filename.c_str(), 0/*world rank*/, FLOAT, &get_rank, &get_dims, &num_chunks, &get_chunk_rank, &get_chunk_dims);
	
			std::cout << "Rank: " << get_rank << ": ";
			for (size_t i = 0; i < get_rank; i++) std::cout << get_dims[i] << " ";
			std::cout << std::endl;
	
			std::cout << "Num Chunks: " << num_chunks << " Chunk rank: " << get_chunk_rank << " ";
			for (size_t i = 0; i < get_chunk_rank; i++) std::cout << get_chunk_dims[i] << " ";
			std::cout << std::endl;
	
			filename = "data/" + filename;
			size_t total_size = 1;
			for (size_t i = 0; i < get_chunk_rank; i++) total_size *= get_chunk_dims[i];
	
			FILE *fp = fopen(filename.c_str(), "rb");
			buffer = (float*)calloc(sizeof(float), total_size);
			fread(buffer, sizeof(float), total_size, fp);
	
			for (size_t i = 0; i < total_size; i++) {
				if (std::abs(buffer[i] - data[i]) > EPS) {
					std::cerr << i << " " << buffer[i] << " vs " << data[i] << std::endl;
					exit(1);
				}
			}
	
			free(get_dims); free(get_chunk_dims); free(data);
			std::cout << "Verified " << filename.c_str() << std::endl;
		}
	}
}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	int rc;
	load_config();
	rc = aoi_init(clovis_local_addr, clovis_ha_addr, clovis_prof,
			clovis_proc_fid, clovis_block_size, tier);
	assert(rc == 0);
	std::cout << "Connected to mero..." << std::endl;
	//verify_particles(4, 13824, 5);
	verify_quantities(5);
	aoi_fini();
	MPI_Finalize();
	return 0;
}
