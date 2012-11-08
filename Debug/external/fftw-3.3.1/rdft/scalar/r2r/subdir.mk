################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../external/fftw-3.3.1/rdft/scalar/r2r/codlist.o \
../external/fftw-3.3.1/rdft/scalar/r2r/e01_8.o \
../external/fftw-3.3.1/rdft/scalar/r2r/e10_8.o 

C_SRCS += \
../external/fftw-3.3.1/rdft/scalar/r2r/codlist.c \
../external/fftw-3.3.1/rdft/scalar/r2r/e01_8.c \
../external/fftw-3.3.1/rdft/scalar/r2r/e10_8.c 

OBJS += \
./external/fftw-3.3.1/rdft/scalar/r2r/codlist.o \
./external/fftw-3.3.1/rdft/scalar/r2r/e01_8.o \
./external/fftw-3.3.1/rdft/scalar/r2r/e10_8.o 

C_DEPS += \
./external/fftw-3.3.1/rdft/scalar/r2r/codlist.d \
./external/fftw-3.3.1/rdft/scalar/r2r/e01_8.d \
./external/fftw-3.3.1/rdft/scalar/r2r/e10_8.d 


# Each subdirectory must supply rules for building sources it contributes
external/fftw-3.3.1/rdft/scalar/r2r/%.o: ../external/fftw-3.3.1/rdft/scalar/r2r/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


