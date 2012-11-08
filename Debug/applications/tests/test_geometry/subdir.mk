################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../applications/tests/test_geometry/test_geometry_main.o 

CPP_SRCS += \
../applications/tests/test_geometry/test_geometry_main.cpp 

OBJS += \
./applications/tests/test_geometry/test_geometry_main.o 

CPP_DEPS += \
./applications/tests/test_geometry/test_geometry_main.d 


# Each subdirectory must supply rules for building sources it contributes
applications/tests/test_geometry/%.o: ../applications/tests/test_geometry/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


