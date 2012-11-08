################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../applications/programs/template_mpi/template_mpi_main.cpp 

OBJS += \
./applications/programs/template_mpi/template_mpi_main.o 

CPP_DEPS += \
./applications/programs/template_mpi/template_mpi_main.d 


# Each subdirectory must supply rules for building sources it contributes
applications/programs/template_mpi/%.o: ../applications/programs/template_mpi/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


