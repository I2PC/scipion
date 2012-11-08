################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/mpi4py-1.2.2/demo/helloworld.c 

CXX_SRCS += \
../external/python/mpi4py-1.2.2/demo/helloworld.cxx 

OBJS += \
./external/python/mpi4py-1.2.2/demo/helloworld.o 

C_DEPS += \
./external/python/mpi4py-1.2.2/demo/helloworld.d 

CXX_DEPS += \
./external/python/mpi4py-1.2.2/demo/helloworld.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/mpi4py-1.2.2/demo/%.o: ../external/python/mpi4py-1.2.2/demo/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/python/mpi4py-1.2.2/demo/%.o: ../external/python/mpi4py-1.2.2/demo/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


