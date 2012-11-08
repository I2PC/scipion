################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../external/fftw-3.3.1/dft/scalar/n.o \
../external/fftw-3.3.1/dft/scalar/t.o 

C_SRCS += \
../external/fftw-3.3.1/dft/scalar/n.c \
../external/fftw-3.3.1/dft/scalar/t.c 

OBJS += \
./external/fftw-3.3.1/dft/scalar/n.o \
./external/fftw-3.3.1/dft/scalar/t.o 

C_DEPS += \
./external/fftw-3.3.1/dft/scalar/n.d \
./external/fftw-3.3.1/dft/scalar/t.d 


# Each subdirectory must supply rules for building sources it contributes
external/fftw-3.3.1/dft/scalar/%.o: ../external/fftw-3.3.1/dft/scalar/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


