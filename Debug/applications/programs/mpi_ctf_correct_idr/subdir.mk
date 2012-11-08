################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../applications/programs/mpi_ctf_correct_idr/mpi_ctf_correct_idr_main.o 

CPP_SRCS += \
../applications/programs/mpi_ctf_correct_idr/mpi_ctf_correct_idr_main.cpp 

OBJS += \
./applications/programs/mpi_ctf_correct_idr/mpi_ctf_correct_idr_main.o 

CPP_DEPS += \
./applications/programs/mpi_ctf_correct_idr/mpi_ctf_correct_idr_main.d 


# Each subdirectory must supply rules for building sources it contributes
applications/programs/mpi_ctf_correct_idr/%.o: ../applications/programs/mpi_ctf_correct_idr/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


