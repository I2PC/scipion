################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../external/fftw-3.3.1/tests/bench-bench.o \
../external/fftw-3.3.1/tests/bench-fftw-bench.o \
../external/fftw-3.3.1/tests/bench-hook.o 

C_SRCS += \
../external/fftw-3.3.1/tests/bench.c \
../external/fftw-3.3.1/tests/fftw-bench.c \
../external/fftw-3.3.1/tests/hook.c 

OBJS += \
./external/fftw-3.3.1/tests/bench.o \
./external/fftw-3.3.1/tests/fftw-bench.o \
./external/fftw-3.3.1/tests/hook.o 

C_DEPS += \
./external/fftw-3.3.1/tests/bench.d \
./external/fftw-3.3.1/tests/fftw-bench.d \
./external/fftw-3.3.1/tests/hook.d 


# Each subdirectory must supply rules for building sources it contributes
external/fftw-3.3.1/tests/%.o: ../external/fftw-3.3.1/tests/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


