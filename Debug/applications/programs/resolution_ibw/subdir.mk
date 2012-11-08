################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../applications/programs/resolution_ibw/resolution_ibw_main.cpp 

OBJS += \
./applications/programs/resolution_ibw/resolution_ibw_main.o 

CPP_DEPS += \
./applications/programs/resolution_ibw/resolution_ibw_main.d 


# Each subdirectory must supply rules for building sources it contributes
applications/programs/resolution_ibw/%.o: ../applications/programs/resolution_ibw/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


