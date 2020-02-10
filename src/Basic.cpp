#include "Basic.h"
#include <math.h>

/** method to calculate the parallel dot product with vect1, vect2 having the
 * ghost cells*/
double dotP(FPfield *vect1, FPfield *vect2, int n) {
  double result = 0;
  double local_result = 0;
  for (int i = 0; i < n; i++) local_result += vect1[i] * vect2[i];
  result = local_result;
  return (result);
}
/** method to calculate dot product */
double dot(FPfield *vect1, FPfield *vect2, int n) {
  double result = 0;
  for (int i = 0; i < n; i++) result += vect1[i] * vect2[i];
  return (result);
}
/** method to calculate the square norm of a vector */
double norm2(FPfield **vect, int nx, int ny) {
  double result = 0;
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++) result += vect[i][j] * vect[i][j];
  return (result);
}
/** method to calculate the square norm of a vector */
double norm2(FPfield ***vect, int nx, int ny) {
  double result = 0;
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++) result += vect[i][j][0] * vect[i][j][0];
  return (result);
}

/** method to calculate the square norm of a vector */
double norm2(FPfield *vect, int nx) {
  double result = 0;
  for (int i = 0; i < nx; i++) result += vect[i] * vect[i];
  return (result);
}

/** method to calculate the parallel dot product */
double norm2P(FPfield ***vect, int nx, int ny, int nz) {
  double result = 0;
  double local_result = 0;
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        local_result += vect[i][j][k] * vect[i][j][k];

  // MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM,
  // MPI_COMM_WORLD);
  result = local_result;
  return (result);
}
/** method to calculate the parallel norm of a vector on different processors
 * with the ghost cell */
double norm2P(FPfield *vect, int n) {
  double result = 0;
  double local_result = 0;
  for (int i = 0; i < n; i++) local_result += vect[i] * vect[i];
  // MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM,
  // MPI_COMM_WORLD);
  result = local_result;
  return (result);
}
/** method to calculate the parallel norm of a vector on different processors
 * with the gost cell*/
double normP(FPfield *vect, int n) {
  FPfield result = 0.0;
  FPfield local_result = 0.0;
  for (int i = 0; i < n; i++) local_result += vect[i] * vect[i];

  // MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM,
  // MPI_COMM_WORLD);
  result = local_result;
  return (sqrt(result));
}
/** method to calculate the difference of two vectors*/
void sub(FPfield *res, FPfield *vect1, FPfield *vect2, int n) {
  for (int i = 0; i < n; i++) res[i] = vect1[i] - vect2[i];
}
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
void sum(FPfield *vect1, FPfield *vect2, int n) {
  for (int i = 0; i < n; i++) vect1[i] += vect2[i];
}
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
void sum(FPfield ***vect1, FPfield ***vect2, int nx, int ny, int nz) {
#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++) vect1[i][j][k] += vect2[i][j][k];
}

/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
void sum(FPfield ***vect1, FPfield ****vect2, int nx, int ny, int nz, int ns) {
#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++) vect1[i][j][k] += vect2[ns][i][j][k];
}

/** method to calculate the subtraction of two vectors vector1 = vector1 -
 * vector2*/
void sub(FPfield ***vect1, FPfield ***vect2, int nx, int ny, int nz) {
#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++) vect1[i][j][k] -= vect2[i][j][k];
}

/** method to sum 4 vectors vector1 = alfa*vector1 + beta*vector2 +
 * gamma*vector3 + delta*vector4 */
void sum4(double ***vect1, double alfa, double ***vect2, double beta,
          double ***vect3, double gamma, double ***vect4, double delta,
          double ***vect5, int nx, int ny, int nz) {
#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] =
            alfa * (vect2[i][j][k] + beta * vect3[i][j][k] +
                    gamma * vect4[i][j][k] + delta * vect5[i][j][k]);
}
/** method to calculate the scalar-vector product */
void scale(FPfield *vect, FPfield alfa, int n) {
#pragma omp parallel for
  for (int i = 0; i < n; i++) vect[i] *= alfa;
}

/** method to calculate the scalar-vector product */
void scale(FPfield ***vect, FPfield alfa, int nx, int ny, int nz) {
#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++) vect[i][j][k] *= alfa;
}
/** method to calculate the scalar product */
void scale(double vect[][2][2], double alfa, int nx, int ny, int nz) {
#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++) vect[i][j][k] *= alfa;
}
/** method to calculate the scalar-vector product */
void scale(FPfield ***vect1, FPfield ***vect2, FPfield alfa, int nx, int ny,
           int nz) {
#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++) vect1[i][j][k] = vect2[i][j][k] * alfa;
}

/** method to calculate the scalar-vector product */
void scale(FPfield *vect1, FPfield *vect2, double alfa, int n) {
#pragma omp parallel for
  for (int i = 0; i < n; i++) vect1[i] = vect2[i] * alfa;
}

/** method to calculate vector1 = vector1 + alfa*vector2   */
void addscale(FPfield alfa, FPfield ***vect1, FPfield ***vect2, int nx, int ny,
              int nz) {
#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = vect1[i][j][k] + alfa * vect2[i][j][k];
}
/** add scale for weights */
void addscale(double alfa, double vect1[][2][2], double vect2[][2][2], int nx,
              int ny, int nz) {
#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = vect1[i][j][k] + alfa * vect2[i][j][k];
}

/** method to calculate vector1 = vector1 + alfa*vector2   */
void addscale(FPfield alfa, FPfield *vect1, FPfield *vect2, int n) {
#pragma omp parallel for
  for (int i = 0; i < n; i++) vect1[i] += alfa * vect2[i];
}
/** method to calculate vector1 = beta*vector1 + alfa*vector2   */
void addscale(FPfield alfa, FPfield beta, FPfield *vect1, FPfield *vect2,
              int n) {
#pragma omp parallel for
  for (int i = 0; i < n; i++) vect1[i] = vect1[i] * beta + alfa * vect2[i];
}
/** method to calculate vector1 = beta*vector1 + alfa*vector2 */
void addscale(FPfield alfa, FPfield beta, FPfield ***vect1, FPfield ***vect2,
              int nx, int ny, int nz) {
#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++) {
          vect1[i][j][k] = beta * vect1[i][j][k] + alfa * vect2[i][j][k];
        }
}

/** method to calculate vector1 = alfa*vector2 + beta*vector3 */
void scaleandsum(double ***vect1, double alfa, double beta, double ***vect2,
                 double ***vect3, int nx, int ny, int nz) {
#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = alfa * vect2[i][j][k] + beta * vect3[i][j][k];
}
/** method to calculate vector1 = alfa*vector2 + beta*vector3 with vector2
 * depending on species*/
void scaleandsum(double ***vect1, double alfa, double beta, double ****vect2,
                 double ***vect3, int ns, int nx, int ny, int nz) {
#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = alfa * vect2[ns][i][j][k] + beta * vect3[i][j][k];
}
/** method to calculate vector1 = alfa*vector2*vector3 with vector2 depending on
 * species*/
void prod(double ***vect1, double alfa, double ****vect2, int ns,
          double ***vect3, int nx, int ny, int nz) {
#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = alfa * vect2[ns][i][j][k] * vect3[i][j][k];
}
/** method to calculate vect1 = vect2/alfa */
void div(double ***vect1, double alfa, double ***vect2, int nx, int ny,
         int nz) {
#pragma omp parallel for
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++) vect1[i][j][k] = vect2[i][j][k] / alfa;
}
void prod6(double ***vect1, double ***vect2, double ***vect3, double ***vect4,
           double ***vect5, double ***vect6, double ***vect7, int nx, int ny,
           int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = vect2[i][j][k] * vect3[i][j][k] +
                         vect4[i][j][k] * vect5[i][j][k] +
                         vect6[i][j][k] * vect7[i][j][k];
}
/** method used for calculating PI */
void proddiv(double ***vect1, double ***vect2, double alfa, double ***vect3,
             double ***vect4, double ***vect5, double ***vect6, double beta,
             double ***vect7, double ***vect8, double gamma, double ***vect9,
             int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = (vect2[i][j][k] +
                          alfa * (vect3[i][j][k] * vect4[i][j][k] -
                                  vect5[i][j][k] * vect6[i][j][k]) +
                          beta * vect7[i][j][k] * vect8[i][j][k]) /
                         (1 + gamma * vect9[i][j][k]);

  // questo mi convince veramente poco!!!!!!!!!!!!!! CAZZO!!!!!!!!!!!!!!!!!!
  // ***vect1++ = (***vect2++ + alfa*((***vect3++)*(***vect4++) -
  // (***vect5++)*(***vect6++)) +
  // beta*(***vect7++)*(***vect8++))/(1+gamma*(***vect9++));
}
/** method to calculate the opposite of a vector */
void neg(FPfield ***vect, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++) vect[i][j][k] = -vect[i][j][k];
}

/** method to calculate the opposite of a vector */
void neg(FPfield *vect, int n) {
  for (int i = 0; i < n; i++) vect[i] = -vect[i];
}
/** method to set equal two vectors */
void eq(double ***vect1, double ***vect2, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++) vect1[i][j][k] = vect2[i][j][k];
}

/** method to set equal two vectors */
void eq(double ****vect1, double ***vect2, int nx, int ny, int nz, int is) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++) vect1[is][i][j][k] = vect2[i][j][k];
}

void eq(double *vect1, double *vect2, int n) {
  for (int i = 0; i < n; i++) vect1[i] = vect2[i];
}
/** method to set a vector to a Value */
void eqValue(FPfield value, FPfield ***vect, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++) vect[i][j][k] = value;
}
void eqValue(double value, double vect[][2][2], int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
#pragma clang loop vectorize(enable)
      for (int k = 0; k < nz; k++) vect[i][j][k] = value;
}

/** method to set a vector to a Value */
void eqValue(FPfield value, FPfield *vect, int n) {
  for (int i = 0; i < n; i++) vect[i] = value;
}
/** method to put a column in a matrix 2D */
void putColumn(FPfield **Matrix, FPfield *vect, int column, int n) {
  for (int i = 0; i < n; i++) Matrix[i][column] = vect[i];
}
/** method to get a column in a matrix 2D */
void getColumn(FPfield *vect, FPfield **Matrix, int column, int n) {
  for (int i = 0; i < n; i++) vect[i] = Matrix[i][column];
}
/** RIFAI QUESTA PARTE questo e' la tomba delle performance*/
void MODULO(double *x, double L) { *x = *x - floor(*x / L) * L; }
/** method to calculate the epsilon machine */
double eps() {
  double eps;
  int i = 1;
  double num = 1;
  double newsum = 1;
  double oldsum = 1;
  while (true) {
      num = num / (2 * i);
      newsum += num;
      if (newsum == oldsum) break;
      oldsum = newsum;
      i++;
    }
  eps = num * 2;
  return (eps);
}

/** method to convert a 1D field in a 3D field not considering guard cells*/
void solver2phys(FPfield ***vectPhys, FPfield *vectSolver, int nx, int ny,
                 int nz) {
  for (int i = 1; i < nx - 1; i++)
    for (int j = 1; j < ny - 1; j++)
      for (int k = 1; k < nz - 1; k++) vectPhys[i][j][k] = *vectSolver++;
}
/** method to convert a 1D field in a 3D field not considering guard cells*/
void solver2phys(FPfield ***vectPhys1, FPfield ***vectPhys2,
                 FPfield ***vectPhys3, FPfield *vectSolver, int nx, int ny,
                 int nz) {
  for (int i = 1; i < nx - 1; i++)
    for (int j = 1; j < ny - 1; j++)
      for (int k = 1; k < nz - 1; k++) {
        vectPhys1[i][j][k] = *vectSolver++;
        vectPhys2[i][j][k] = *vectSolver++;
        vectPhys3[i][j][k] = *vectSolver++;
      }
}

/** method to convert a 3D field in a 1D field not considering guard cells*/
void phys2solver(FPfield *vectSolver, FPfield ***vectPhys, int nx, int ny,
                 int nz) {
  for (int i = 1; i < nx - 1; i++)
    for (int j = 1; j < ny - 1; j++)
      for (int k = 1; k < nz - 1; k++) *vectSolver++ = vectPhys[i][j][k];
}
/** method to convert a 3D field in a 1D field not considering guard cells*/
void phys2solver(FPfield *vectSolver, FPfield ***vectPhys1,
                 FPfield ***vectPhys2, FPfield ***vectPhys3, int nx, int ny,
                 int nz) {
  for (int i = 1; i < nx - 1; i++)
    for (int j = 1; j < ny - 1; j++)
      for (int k = 1; k < nz - 1; k++) {
        *vectSolver++ = vectPhys1[i][j][k];
        *vectSolver++ = vectPhys2[i][j][k];
        *vectSolver++ = vectPhys3[i][j][k];
      }
}
