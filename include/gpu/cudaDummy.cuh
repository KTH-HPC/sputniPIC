// CREDIT TO: https://stackoverflow.com/questions/51959774/cuda-in-vscode

#pragma once
#ifdef __INTELLISENSE__
void __syncthreads();  // workaround __syncthreads warning
#define KERNEL_ARG2(grid, block)
#define KERNEL_ARG3(grid, block, sh_mem)
#define KERNEL_ARG4(grid, block, sh_mem, stream)
#else
#define KERNEL_ARG2(grid, block) <<< grid, block >>>
#define KERNEL_ARG3(grid, block, sh_mem) <<< grid, block, sh_mem >>>
#define KERNEL_ARG4(grid, block, sh_mem, stream) \
<<< grid, block, sh_mem, stream >>>
#endif