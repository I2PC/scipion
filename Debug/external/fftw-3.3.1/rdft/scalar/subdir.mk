################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../external/fftw-3.3.1/rdft/scalar/hc2c.o \
../external/fftw-3.3.1/rdft/scalar/hfb.o \
../external/fftw-3.3.1/rdft/scalar/r2c.o \
../external/fftw-3.3.1/rdft/scalar/r2r.o 

C_SRCS += \
../external/fftw-3.3.1/rdft/scalar/hc2c.c \
../external/fftw-3.3.1/rdft/scalar/hfb.c \
../external/fftw-3.3.1/rdft/scalar/r2c.c \
../external/fftw-3.3.1/rdft/scalar/r2r.c 

OBJS += \
./external/fftw-3.3.1/rdft/scalar/hc2c.o \
./external/fftw-3.3.1/rdft/scalar/hfb.o \
./external/fftw-3.3.1/rdft/scalar/r2c.o \
./external/fftw-3.3.1/rdft/scalar/r2r.o 

C_DEPS += \
./external/fftw-3.3.1/rdft/scalar/hc2c.d \
./external/fftw-3.3.1/rdft/scalar/hfb.d \
./external/fftw-3.3.1/rdft/scalar/r2c.d \
./external/fftw-3.3.1/rdft/scalar/r2r.d 


# Each subdirectory must supply rules for building sources it contributes
external/fftw-3.3.1/rdft/scalar/%.o: ../external/fftw-3.3.1/rdft/scalar/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


