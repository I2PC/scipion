################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../applications/programs/angular_continuous_assign/angular_continuous_assign_main.o 

CPP_SRCS += \
../applications/programs/angular_continuous_assign/angular_continuous_assign_main.cpp 

OBJS += \
./applications/programs/angular_continuous_assign/angular_continuous_assign_main.o 

CPP_DEPS += \
./applications/programs/angular_continuous_assign/angular_continuous_assign_main.d 


# Each subdirectory must supply rules for building sources it contributes
applications/programs/angular_continuous_assign/%.o: ../applications/programs/angular_continuous_assign/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


