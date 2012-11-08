################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/numpy-1.6.1/numpy/linalg/blas_lite.c \
../external/python/numpy-1.6.1/numpy/linalg/dlamch.c \
../external/python/numpy-1.6.1/numpy/linalg/dlapack_lite.c \
../external/python/numpy-1.6.1/numpy/linalg/f2c_lite.c \
../external/python/numpy-1.6.1/numpy/linalg/lapack_litemodule.c \
../external/python/numpy-1.6.1/numpy/linalg/python_xerbla.c \
../external/python/numpy-1.6.1/numpy/linalg/zlapack_lite.c 

OBJS += \
./external/python/numpy-1.6.1/numpy/linalg/blas_lite.o \
./external/python/numpy-1.6.1/numpy/linalg/dlamch.o \
./external/python/numpy-1.6.1/numpy/linalg/dlapack_lite.o \
./external/python/numpy-1.6.1/numpy/linalg/f2c_lite.o \
./external/python/numpy-1.6.1/numpy/linalg/lapack_litemodule.o \
./external/python/numpy-1.6.1/numpy/linalg/python_xerbla.o \
./external/python/numpy-1.6.1/numpy/linalg/zlapack_lite.o 

C_DEPS += \
./external/python/numpy-1.6.1/numpy/linalg/blas_lite.d \
./external/python/numpy-1.6.1/numpy/linalg/dlamch.d \
./external/python/numpy-1.6.1/numpy/linalg/dlapack_lite.d \
./external/python/numpy-1.6.1/numpy/linalg/f2c_lite.d \
./external/python/numpy-1.6.1/numpy/linalg/lapack_litemodule.d \
./external/python/numpy-1.6.1/numpy/linalg/python_xerbla.d \
./external/python/numpy-1.6.1/numpy/linalg/zlapack_lite.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/numpy-1.6.1/numpy/linalg/%.o: ../external/python/numpy-1.6.1/numpy/linalg/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


