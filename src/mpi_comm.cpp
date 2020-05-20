#include "mpi_comm.h"


template <typename T>
MPI_Datatype _mpi_get_datatype(){
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


void _reduce_copy(FPinterp* array, FPinterp* recv_buf, int length){

	MPI_Datatype type = _mpi_get_datatype<FPinterp>();
	MPI_Reduce(array, recv_buf, length, type, MPI_SUM, 0, MPI_COMM_WORLD);
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0){
		memcpy(array, recv_buf, sizeof(FPinterp)*length);
	}
}

void mpi_reduce_densities(struct grid* grd, struct interpDensNet* idn){

	int count = grd->nxn*grd->nyn*grd->nzn;
	FPinterp* recv_buf = (FPinterp*) malloc(sizeof(FPinterp)*count);

	_reduce_copy(idn->rhon_flat, recv_buf, count);
	_reduce_copy(idn->rhoc_flat, recv_buf, grd->nxc*grd->nyc*grd->nzc);

	_reduce_copy(idn->Jx_flat, recv_buf, count);
	_reduce_copy(idn->Jy_flat, recv_buf, count);
	_reduce_copy(idn->Jz_flat, recv_buf, count);

	_reduce_copy(idn->pxx_flat, recv_buf, count);
	_reduce_copy(idn->pxy_flat, recv_buf, count);
	_reduce_copy(idn->pxz_flat, recv_buf, count);

	_reduce_copy(idn->pyy_flat, recv_buf, count);
	_reduce_copy(idn->pyz_flat, recv_buf, count);
	_reduce_copy(idn->pzz_flat, recv_buf, count);

	free(recv_buf);

}


void mpi_reduce_densities(struct grid* grd, struct interpDensSpecies* ids){

	int count = grd->nxn*grd->nyn*grd->nzn;
	FPinterp* recv_buf = (FPinterp*) malloc(sizeof(FPinterp)*count);

	_reduce_copy(ids->rhon_flat, recv_buf, count);
	_reduce_copy(ids->rhoc_flat, recv_buf, grd->nxc*grd->nyc*grd->nzc);

	_reduce_copy(ids->Jx_flat, recv_buf, count);
	_reduce_copy(ids->Jy_flat, recv_buf, count);
	_reduce_copy(ids->Jz_flat, recv_buf, count);

	_reduce_copy(ids->pxx_flat, recv_buf, count);
	_reduce_copy(ids->pxy_flat, recv_buf, count);
	_reduce_copy(ids->pxz_flat, recv_buf, count);

	_reduce_copy(ids->pyy_flat, recv_buf, count);
	_reduce_copy(ids->pyz_flat, recv_buf, count);
	_reduce_copy(ids->pzz_flat, recv_buf, count);

	free(recv_buf);

}


void mpi_broadcast_field(struct grid *grd, struct EMfield *field){

	int count = grd->nxn*grd->nyn*grd->nzn;

	MPI_Datatype type = _mpi_get_datatype<FPinterp>();

	MPI_Bcast(field->Ex_flat, count, type, 0, MPI_COMM_WORLD);
	MPI_Bcast(field->Ey_flat, count, type, 0, MPI_COMM_WORLD);
	MPI_Bcast(field->Ez_flat, count, type, 0, MPI_COMM_WORLD);

	MPI_Bcast(field->Bxn_flat, count, type, 0, MPI_COMM_WORLD);
	MPI_Bcast(field->Byn_flat, count, type, 0, MPI_COMM_WORLD);
	MPI_Bcast(field->Bzn_flat, count, type, 0, MPI_COMM_WORLD);


}