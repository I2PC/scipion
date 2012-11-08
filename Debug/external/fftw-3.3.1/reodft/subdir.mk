################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../external/fftw-3.3.1/reodft/conf.o \
../external/fftw-3.3.1/reodft/redft00e-r2hc-pad.o \
../external/fftw-3.3.1/reodft/reodft00e-splitradix.o \
../external/fftw-3.3.1/reodft/reodft010e-r2hc.o \
../external/fftw-3.3.1/reodft/reodft11e-r2hc-odd.o \
../external/fftw-3.3.1/reodft/reodft11e-radix2.o \
../external/fftw-3.3.1/reodft/rodft00e-r2hc-pad.o 

C_SRCS += \
../external/fftw-3.3.1/reodft/conf.c \
../external/fftw-3.3.1/reodft/redft00e-r2hc-pad.c \
../external/fftw-3.3.1/reodft/redft00e-r2hc.c \
../external/fftw-3.3.1/reodft/reodft00e-splitradix.c \
../external/fftw-3.3.1/reodft/reodft010e-r2hc.c \
../external/fftw-3.3.1/reodft/reodft11e-r2hc-odd.c \
../external/fftw-3.3.1/reodft/reodft11e-r2hc.c \
../external/fftw-3.3.1/reodft/reodft11e-radix2.c \
../external/fftw-3.3.1/reodft/rodft00e-r2hc-pad.c \
../external/fftw-3.3.1/reodft/rodft00e-r2hc.c 

OBJS += \
./external/fftw-3.3.1/reodft/conf.o \
./external/fftw-3.3.1/reodft/redft00e-r2hc-pad.o \
./external/fftw-3.3.1/reodft/redft00e-r2hc.o \
./external/fftw-3.3.1/reodft/reodft00e-splitradix.o \
./external/fftw-3.3.1/reodft/reodft010e-r2hc.o \
./external/fftw-3.3.1/reodft/reodft11e-r2hc-odd.o \
./external/fftw-3.3.1/reodft/reodft11e-r2hc.o \
./external/fftw-3.3.1/reodft/reodft11e-radix2.o \
./external/fftw-3.3.1/reodft/rodft00e-r2hc-pad.o \
./external/fftw-3.3.1/reodft/rodft00e-r2hc.o 

C_DEPS += \
./external/fftw-3.3.1/reodft/conf.d \
./external/fftw-3.3.1/reodft/redft00e-r2hc-pad.d \
./external/fftw-3.3.1/reodft/redft00e-r2hc.d \
./external/fftw-3.3.1/reodft/reodft00e-splitradix.d \
./external/fftw-3.3.1/reodft/reodft010e-r2hc.d \
./external/fftw-3.3.1/reodft/reodft11e-r2hc-odd.d \
./external/fftw-3.3.1/reodft/reodft11e-r2hc.d \
./external/fftw-3.3.1/reodft/reodft11e-radix2.d \
./external/fftw-3.3.1/reodft/rodft00e-r2hc-pad.d \
./external/fftw-3.3.1/reodft/rodft00e-r2hc.d 


# Each subdirectory must supply rules for building sources it contributes
external/fftw-3.3.1/reodft/%.o: ../external/fftw-3.3.1/reodft/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


