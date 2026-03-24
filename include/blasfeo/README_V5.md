# BLASFEO Build for VEX V5 (ARMv7-A)

This document explains how to resolve the "undefined reference" issues in BLASFEO when targeting the VEX V5.

## The Problem
BLASFEO's optimized assembly kernels are guarded by macros like `OS_LINUX` or `OS_MAC`. Since the PROS toolchain targets `arm-none-eabi` (bare-metal ELF), these macros are not defined by default, causing the kernels to be skipped during compilation. This results in missing math symbols in `libblasfeo.a`.

## The Fix

### 1. Global Flags
Add `-DOS_LINUX` to your `common.mk` or `Makefile` to enable the kernel assembly code:

```makefile
CPPFLAGS += -DOS_LINUX
```

### 2. Manual Kernel Compilation
If `libblasfeo.a` is missing symbols (e.g., `kernel_sgemm_nt_8x4_lib44cc`), you must manually compile the high-performance kernels and integrate them.

**Important:** Use `MACRO_LEVEL=1` to ensure internal kernel references are correctly linked.

#### Compilation Command Template:
```bash
arm-none-eabi-gcc -c -O2 -mcpu=cortex-a9 -marm -mfloat-abi=hard -mfpu=neon \
    -DOS_LINUX -DTARGET_ARMV7A_ARM_CORTEX_A9 -DLA_HIGH_PERFORMANCE -DMF_PANELMAJ \
    -DBLAS_API -DMACRO_LEVEL=1 -Iinclude -Ikernel/armv7a \
    -D"PROLOGUE=stmdb sp!, {r4 - r10, fp, lr}; add fp, sp, #36; fstmfdd sp!, {d8-d15};" \
    -D"EPILOGUE=fldmfdd sp!, {d8-d15}; ldmia sp!, {r4 - r10, fp, pc};" \
    kernel/armv7a/[KERNEL_FILE].S -o [KERNEL_FILE].o
```

### 3. Library Integration
Use `ar` to update `firmware/libblasfeo.a`. 

**Note on Symbol Conflicts:** 
The `lib4.S` kernels (e.g., `kernel_sgemm_4x4_lib4.S`) already include the base `lib.S` code. If you compile both as separate objects, you will get "multiple definition" errors. 

**Correct Workflow:**
1. Compile only the `lib4` variants.
2. Remove any existing base `lib` objects from the archive:
   ```bash
   arm-none-eabi-ar d firmware/libblasfeo.a kernel_sgemm_8x4_lib.o kernel_sgemm_4x4_lib.o kernel_dgemm_4x4_lib.o
   ```
3. Add the new `lib4` objects:
   ```bash
   arm-none-eabi-ar rcs firmware/libblasfeo.a [COMPILED_OBJECTS]
   ```

## Essential Kernels for V5
The following kernels are usually required for high-performance BLAS/LAPACK operations:
- `kernel_dgemm_4x4_lib4.S`
- `kernel_sgemm_8x4_lib4.S`
- `kernel_sgemm_4x4_lib4.S`
- `kernel_sgemm_12x4_lib4.S`
