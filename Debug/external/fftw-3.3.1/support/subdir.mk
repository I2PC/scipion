################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/fftw-3.3.1/support/addchain.c 

OBJS += \
./external/fftw-3.3.1/support/addchain.o 

C_DEPS += \
./external/fftw-3.3.1/support/addchain.d 


# Each subdirectory must supply rules for building sources it contributes
external/fftw-3.3.1/support/%.o: ../external/fftw-3.3.1/support/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


