################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../applications/programs/mrc_create_metadata/mrc_create_metadata_main.o 

CPP_SRCS += \
../applications/programs/mrc_create_metadata/mrc_create_metadata_main.cpp 

OBJS += \
./applications/programs/mrc_create_metadata/mrc_create_metadata_main.o 

CPP_DEPS += \
./applications/programs/mrc_create_metadata/mrc_create_metadata_main.d 


# Each subdirectory must supply rules for building sources it contributes
applications/programs/mrc_create_metadata/%.o: ../applications/programs/mrc_create_metadata/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


