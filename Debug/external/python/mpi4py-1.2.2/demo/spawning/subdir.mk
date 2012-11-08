################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/mpi4py-1.2.2/demo/spawning/cpi-master.c \
../external/python/mpi4py-1.2.2/demo/spawning/cpi-worker.c 

CXX_SRCS += \
../external/python/mpi4py-1.2.2/demo/spawning/cpi-master.cxx \
../external/python/mpi4py-1.2.2/demo/spawning/cpi-worker.cxx 

OBJS += \
./external/python/mpi4py-1.2.2/demo/spawning/cpi-master.o \
./external/python/mpi4py-1.2.2/demo/spawning/cpi-worker.o 

C_DEPS += \
./external/python/mpi4py-1.2.2/demo/spawning/cpi-master.d \
./external/python/mpi4py-1.2.2/demo/spawning/cpi-worker.d 

CXX_DEPS += \
./external/python/mpi4py-1.2.2/demo/spawning/cpi-master.d \
./external/python/mpi4py-1.2.2/demo/spawning/cpi-worker.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/mpi4py-1.2.2/demo/spawning/%.o: ../external/python/mpi4py-1.2.2/demo/spawning/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/python/mpi4py-1.2.2/demo/spawning/%.o: ../external/python/mpi4py-1.2.2/demo/spawning/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


