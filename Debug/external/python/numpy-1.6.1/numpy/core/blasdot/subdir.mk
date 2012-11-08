################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/numpy-1.6.1/numpy/core/blasdot/_dotblas.c 

OBJS += \
./external/python/numpy-1.6.1/numpy/core/blasdot/_dotblas.o 

C_DEPS += \
./external/python/numpy-1.6.1/numpy/core/blasdot/_dotblas.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/numpy-1.6.1/numpy/core/blasdot/%.o: ../external/python/numpy-1.6.1/numpy/core/blasdot/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


