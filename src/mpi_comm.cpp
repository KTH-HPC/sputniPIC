#include "mpi_comm.h"


template <typename T>
MPI_Datatype _mpi_get_basetype(){
	if(typeid(T) == typeid(float)){
		return MPI_FLOAT;
	}
	else if(typeid(T) == typeid(double)){
		return MPI_DOUBLE;
	}
	else{
		printf("Unknown precision type '%s'", typeid(T).name());
		exit(EXIT_FAILURE);
	}
}

inline void _scatter(FPpart *send, FPpart *recv, long count, MPI_Datatype type){

	MPI_Scatter(
		send, count, type,
		recv, count, type,
		0, MPI_COMM_WORLD
		);
}

void mpi_scatter_particles(particles *part_global, particles *part_local){

	MPI_Datatype type = _mpi_get_basetype<FPpart>();

	_scatter(part_global->x, part_local->x, part_local->nop, type);
	_scatter(part_global->y, part_local->y, part_local->nop, type);
	_scatter(part_global->z, part_local->z, part_local->nop, type);
	_scatter(part_global->u, part_local->u, part_local->nop, type);
	_scatter(part_global->v, part_local->v, part_local->nop, type);
	_scatter(part_global->w, part_local->w, part_local->nop, type);

	type = _mpi_get_basetype<FPinterp>();
	MPI_Scatter(
		part_global->q, part_local->nop, type,
		part_local->q, part_local->nop, type,
		0, MPI_COMM_WORLD
	);

	MPI_Scatter(
		part_global->track_particle, part_local->nop, MPI_C_BOOL,
		part_local->track_particle, part_local->nop, MPI_C_BOOL,
		0, MPI_COMM_WORLD
	);
}

void mpi_gather_particles(particles *part_global, particles *part_local){

	

}

inline void _reduce_interp(FPinterp* array, int length, int rank){

	MPI_Datatype type = _mpi_get_basetype<FPinterp>();
	void* recv_buf = nullptr;

	if(!rank){
		MPI_Reduce(MPI_IN_PLACE, array, length, type, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	else{
		MPI_Reduce(array, recv_buf, length, type, MPI_SUM, 0, MPI_COMM_WORLD);
	}
}

void mpi_reduce_densities(struct grid* grd, struct interpDensNet* idn){

	int rank;
	int count = grd->nxn*grd->nyn*grd->nzn;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	_reduce_interp(idn->rhon_flat, count, rank);
	_reduce_interp(idn->rhoc_flat, grd->nxc*grd->nyc*grd->nzc, rank);

	_reduce_interp(idn->Jx_flat, count, rank);
	_reduce_interp(idn->Jy_flat, count, rank);
	_reduce_interp(idn->Jz_flat, count, rank);

	_reduce_interp(idn->pxx_flat, count, rank);
	_reduce_interp(idn->pxy_flat, count, rank);
	_reduce_interp(idn->pxz_flat, count, rank);

	_reduce_interp(idn->pyy_flat, count, rank);
	_reduce_interp(idn->pyz_flat, count, rank);
	_reduce_interp(idn->pzz_flat, count, rank);

}

void mpi_reduce_densities(struct grid* grd, struct interpDensSpecies* ids){

	int rank;
	int count = grd->nxn*grd->nyn*grd->nzn;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	_reduce_interp(ids->rhon_flat, count, rank);
	_reduce_interp(ids->rhoc_flat, grd->nxc*grd->nyc*grd->nzc, rank);

	_reduce_interp(ids->Jx_flat, count, rank);
	_reduce_interp(ids->Jy_flat, count, rank);
	_reduce_interp(ids->Jz_flat, count, rank);

	_reduce_interp(ids->pxx_flat, count, rank);
	_reduce_interp(ids->pxy_flat, count, rank);
	_reduce_interp(ids->pxz_flat, count, rank);

	_reduce_interp(ids->pyy_flat, count, rank);
	_reduce_interp(ids->pyz_flat, count, rank);
	_reduce_interp(ids->pzz_flat, count, rank);


}


void mpi_broadcast_field(struct grid *grd, struct EMfield *field){

	int count = grd->nxn*grd->nyn*grd->nzn;

	MPI_Datatype type = _mpi_get_basetype<FPinterp>();

	MPI_Bcast(field->Ex_flat, count, type, 0, MPI_COMM_WORLD);
	MPI_Bcast(field->Ey_flat, count, type, 0, MPI_COMM_WORLD);
	MPI_Bcast(field->Ez_flat, count, type, 0, MPI_COMM_WORLD);

	MPI_Bcast(field->Bxn_flat, count, type, 0, MPI_COMM_WORLD);
	MPI_Bcast(field->Byn_flat, count, type, 0, MPI_COMM_WORLD);
	MPI_Bcast(field->Bzn_flat, count, type, 0, MPI_COMM_WORLD);


}