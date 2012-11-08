################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/mpi4py-1.2.2/src/MPE.c \
../external/python/mpi4py-1.2.2/src/MPI.c \
../external/python/mpi4py-1.2.2/src/dynload.c \
../external/python/mpi4py-1.2.2/src/mpi4py.MPE.c \
../external/python/mpi4py-1.2.2/src/mpi4py.MPI.c \
../external/python/mpi4py-1.2.2/src/python.c 

OBJS += \
./external/python/mpi4py-1.2.2/src/MPE.o \
./external/python/mpi4py-1.2.2/src/MPI.o \
./external/python/mpi4py-1.2.2/src/dynload.o \
./external/python/mpi4py-1.2.2/src/mpi4py.MPE.o \
./external/python/mpi4py-1.2.2/src/mpi4py.MPI.o \
./external/python/mpi4py-1.2.2/src/python.o 

C_DEPS += \
./external/python/mpi4py-1.2.2/src/MPE.d \
./external/python/mpi4py-1.2.2/src/MPI.d \
./external/python/mpi4py-1.2.2/src/dynload.d \
./external/python/mpi4py-1.2.2/src/mpi4py.MPE.d \
./external/python/mpi4py-1.2.2/src/mpi4py.MPI.d \
./external/python/mpi4py-1.2.2/src/python.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/mpi4py-1.2.2/src/%.o: ../external/python/mpi4py-1.2.2/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


