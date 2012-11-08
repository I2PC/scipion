################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../applications/programs/angular_commonline/angular_commonline_main.cpp 

OBJS += \
./applications/programs/angular_commonline/angular_commonline_main.o 

CPP_DEPS += \
./applications/programs/angular_commonline/angular_commonline_main.d 


# Each subdirectory must supply rules for building sources it contributes
applications/programs/angular_commonline/%.o: ../applications/programs/angular_commonline/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


