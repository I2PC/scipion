################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../external/fftw-3.3.1/tools/fftw_wisdom-fftw-wisdom.o 

C_SRCS += \
../external/fftw-3.3.1/tools/fftw-wisdom.c 

OBJS += \
./external/fftw-3.3.1/tools/fftw-wisdom.o 

C_DEPS += \
./external/fftw-3.3.1/tools/fftw-wisdom.d 


# Each subdirectory must supply rules for building sources it contributes
external/fftw-3.3.1/tools/%.o: ../external/fftw-3.3.1/tools/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


