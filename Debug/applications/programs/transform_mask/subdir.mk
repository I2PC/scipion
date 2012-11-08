################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../applications/programs/transform_mask/transform_mask_main.o 

CPP_SRCS += \
../applications/programs/transform_mask/transform_mask_main.cpp 

OBJS += \
./applications/programs/transform_mask/transform_mask_main.o 

CPP_DEPS += \
./applications/programs/transform_mask/transform_mask_main.d 


# Each subdirectory must supply rules for building sources it contributes
applications/programs/transform_mask/%.o: ../applications/programs/transform_mask/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


