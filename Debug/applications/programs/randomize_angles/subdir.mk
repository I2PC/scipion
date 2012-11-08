################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../applications/programs/randomize_angles/randomize_angles_main.cpp 

OBJS += \
./applications/programs/randomize_angles/randomize_angles_main.o 

CPP_DEPS += \
./applications/programs/randomize_angles/randomize_angles_main.d 


# Each subdirectory must supply rules for building sources it contributes
applications/programs/randomize_angles/%.o: ../applications/programs/randomize_angles/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


