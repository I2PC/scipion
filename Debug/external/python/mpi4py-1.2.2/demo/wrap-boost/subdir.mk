################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CXX_SRCS += \
../external/python/mpi4py-1.2.2/demo/wrap-boost/helloworld.cxx 

OBJS += \
./external/python/mpi4py-1.2.2/demo/wrap-boost/helloworld.o 

CXX_DEPS += \
./external/python/mpi4py-1.2.2/demo/wrap-boost/helloworld.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/mpi4py-1.2.2/demo/wrap-boost/%.o: ../external/python/mpi4py-1.2.2/demo/wrap-boost/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


