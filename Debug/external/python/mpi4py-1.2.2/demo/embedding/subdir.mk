################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/mpi4py-1.2.2/demo/embedding/helloworld.c 

OBJS += \
./external/python/mpi4py-1.2.2/demo/embedding/helloworld.o 

C_DEPS += \
./external/python/mpi4py-1.2.2/demo/embedding/helloworld.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/mpi4py-1.2.2/demo/embedding/%.o: ../external/python/mpi4py-1.2.2/demo/embedding/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


