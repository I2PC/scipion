################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../applications/programs/mpi_mlf_tomo/mpi_mlf_tomo_main.cpp 

OBJS += \
./applications/programs/mpi_mlf_tomo/mpi_mlf_tomo_main.o 

CPP_DEPS += \
./applications/programs/mpi_mlf_tomo/mpi_mlf_tomo_main.d 


# Each subdirectory must supply rules for building sources it contributes
applications/programs/mpi_mlf_tomo/%.o: ../applications/programs/mpi_mlf_tomo/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


