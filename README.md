
# sputniPIC
sputniPIC is a space plasma simulation software that uses the Particle-In-Cell (PIC) method.
## Building
sputniPIC depends on MPI, OpenMP, and the optionally CUDA.
```bash
$ git clone https://github.com/steven-chien/sputniPIC.git
$ cd sputniPIC
$ make VERSION=GPU
```
If you want to use the CPU version, specify version to CPU instead.
```bash
$ make VERSION=CPU
```
## Usage
```bash
$ make clean
$ make -j
$ mkdir data_3D
$ ./bin/sputniPIC.out ./inputfiles/GEM_3D.inp |& tee gpu_results.txt

-------------------------
sputniPIC Sim. Parameters
-------------------------
Number of species    = 4
Number of particles of species 0 = 65536000	 (MAX = 65536000)  QOM = -64
Number of particles of species 1 = 65536000	 (MAX = 65536000)  QOM = 1
Number of particles of species 2 = 65536000	 (MAX = 65536000)  QOM = -64
Number of particles of species 3 = 65536000	 (MAX = 65536000)  QOM = 1
x-Length                 = 20
y-Length                 = 10
z-Length                 = 10
Number of cells (x)      = 128
Number of cells (y)      = 64
Number of cells (z)      = 64
Time step                = 0.25
Number of cycles         = 100
Results saved in: data_3D
Output directory data_3D exists.
*************************************************
**  Initialize GEM Challenge with Pertubation  **
*************************************************
** B0x = 0.0195
** B0y = 0
** B0z = 0
** Delta (current sheet thickness) = 0.5
** rho species 0 = 1 CURRENT SHEET 
** rho species 1 = 1 CURRENT SHEET 
** rho species 2 = 0.1 BACKGROUND 
** rho species 3 = 0.1 BACKGROUND 
*************************************************
Writing ic_data...
GPU has free memory: 3941072896x1
CUDA return check: Enabled
Total number of cores: 16
Total number of GPUs: 1
Total number of particles: 262144000; 6000 MB of data.
Allocating 375 MB of memory for particles on gpu
batch_size per species of 4096000 (93 MB)

***********************
   cycle = 1
***********************
***  MOVER  ITERATIONS = 3 - Species 0 *** on gpu 0 Species  - 16 batches 
***  MOVER  ITERATIONS = 3 - Species 1 *** on gpu 0 Species  - 16 batches 
***  MOVER  ITERATIONS = 3 - Species 2 *** on gpu 0 Species  - 16 batches 
***  MOVER  ITERATIONS = 3 - Species 3 *** on gpu 0 Species  - 16 batches 
***********************
*** DIVERGENCE CLEANING ***
Initial error: 0.492688
CG converged at iteration # 72
*** MAXWELL SOLVER ***
Initial residual: 0.546317 norm b vector (source) = 0.823927
GMRES converged at restart # 0; iteration #15 with error: 0.000895017
*** B CALCULATION ***
Timing Cycle 1 : 5.20759 0.0 3.15303 0

***********************
   cycle = 2
***********************
***  MOVER  ITERATIONS = 3 - Species 0 *** on gpu 0 Species  - 16 batches 
***  MOVER  ITERATIONS = 3 - Species 1 *** on gpu 0 Species  - 16 batches 
***  MOVER  ITERATIONS = 3 - Species 2 *** on gpu 0 Species  - 16 batches 
***  MOVER  ITERATIONS = 3 - Species 3 *** on gpu 0 Species  - 16 batches 
***********************
*** DIVERGENCE CLEANING ***
Initial error: 0.519341
CG converged at iteration # 48
*** MAXWELL SOLVER ***
Initial residual: 0.808559 norm b vector (source) = 0.694725
GMRES converged at restart # 0; iteration #13 with error: 0.000898135
*** B CALCULATION ***
Timing Cycle 2 : 4.94224 0.0 2.5049 0


```
## Output format
- Every Iterations
```Timing Cycle <Iterations> : <Mover Time> <Interp Time> <Field Time> <IO Time>```
- Finally
```
Mover: <average> <standard deviation>
Field: <average> <standard deviation>
IO: <average> <standard deviation>
```
## License
The software is released under BSD 2-Clause license. See LICENSE for details.

## Cite us
If you find our implementation and [paper](https://arxiv.org/pdf/2008.04397.pdf) useful to your work, we would apprecite if you cite us with:
```
@INPROCEEDINGS{9235052,
  author={S. W. D. {Chien} and J. {Nylund} and G. {Bengtsson} and I. B. {Peng} and A. {Podobas} and S. {Markidis}},
  booktitle={2020 IEEE 32nd International Symposium on Computer Architecture and High Performance Computing (SBAC-PAD)}, 
  title={sputniPIC: An Implicit Particle-in-Cell Code for Multi-GPU Systems}, 
  year={2020},
  volume={},
  number={},
  pages={149-156},
  doi={10.1109/SBAC-PAD49847.2020.00030}}
```
