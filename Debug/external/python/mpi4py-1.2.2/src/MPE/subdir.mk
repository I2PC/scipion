################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/mpi4py-1.2.2/src/MPE/mpe-log.c \
../external/python/mpi4py-1.2.2/src/MPE/pmpi-ampe.c \
../external/python/mpi4py-1.2.2/src/MPE/pmpi-lmpe.c \
../external/python/mpi4py-1.2.2/src/MPE/pmpi-tmpe.c 

OBJS += \
./external/python/mpi4py-1.2.2/src/MPE/mpe-log.o \
./external/python/mpi4py-1.2.2/src/MPE/pmpi-ampe.o \
./external/python/mpi4py-1.2.2/src/MPE/pmpi-lmpe.o \
./external/python/mpi4py-1.2.2/src/MPE/pmpi-tmpe.o 

C_DEPS += \
./external/python/mpi4py-1.2.2/src/MPE/mpe-log.d \
./external/python/mpi4py-1.2.2/src/MPE/pmpi-ampe.d \
./external/python/mpi4py-1.2.2/src/MPE/pmpi-lmpe.d \
./external/python/mpi4py-1.2.2/src/MPE/pmpi-tmpe.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/mpi4py-1.2.2/src/MPE/%.o: ../external/python/mpi4py-1.2.2/src/MPE/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


