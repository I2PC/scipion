################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../applications/programs/mpi_tomo_extract_subvolume/mpi_tomo_extract_subvolume.o 

CPP_SRCS += \
../applications/programs/mpi_tomo_extract_subvolume/mpi_tomo_extract_subvolume.cpp 

OBJS += \
./applications/programs/mpi_tomo_extract_subvolume/mpi_tomo_extract_subvolume.o 

CPP_DEPS += \
./applications/programs/mpi_tomo_extract_subvolume/mpi_tomo_extract_subvolume.d 


# Each subdirectory must supply rules for building sources it contributes
applications/programs/mpi_tomo_extract_subvolume/%.o: ../applications/programs/mpi_tomo_extract_subvolume/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


