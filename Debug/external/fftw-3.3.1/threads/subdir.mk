################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../external/fftw-3.3.1/threads/libfftw3_threads_la-api.o \
../external/fftw-3.3.1/threads/libfftw3_threads_la-conf.o \
../external/fftw-3.3.1/threads/libfftw3_threads_la-ct.o \
../external/fftw-3.3.1/threads/libfftw3_threads_la-dft-vrank-geq1.o \
../external/fftw-3.3.1/threads/libfftw3_threads_la-f77api.o \
../external/fftw-3.3.1/threads/libfftw3_threads_la-hc2hc.o \
../external/fftw-3.3.1/threads/libfftw3_threads_la-rdft-vrank-geq1.o \
../external/fftw-3.3.1/threads/libfftw3_threads_la-threads.o \
../external/fftw-3.3.1/threads/libfftw3_threads_la-vrank-geq1-rdft2.o 

C_SRCS += \
../external/fftw-3.3.1/threads/api.c \
../external/fftw-3.3.1/threads/conf.c \
../external/fftw-3.3.1/threads/ct.c \
../external/fftw-3.3.1/threads/dft-vrank-geq1.c \
../external/fftw-3.3.1/threads/f77api.c \
../external/fftw-3.3.1/threads/hc2hc.c \
../external/fftw-3.3.1/threads/openmp.c \
../external/fftw-3.3.1/threads/rdft-vrank-geq1.c \
../external/fftw-3.3.1/threads/threads.c \
../external/fftw-3.3.1/threads/vrank-geq1-rdft2.c 

OBJS += \
./external/fftw-3.3.1/threads/api.o \
./external/fftw-3.3.1/threads/conf.o \
./external/fftw-3.3.1/threads/ct.o \
./external/fftw-3.3.1/threads/dft-vrank-geq1.o \
./external/fftw-3.3.1/threads/f77api.o \
./external/fftw-3.3.1/threads/hc2hc.o \
./external/fftw-3.3.1/threads/openmp.o \
./external/fftw-3.3.1/threads/rdft-vrank-geq1.o \
./external/fftw-3.3.1/threads/threads.o \
./external/fftw-3.3.1/threads/vrank-geq1-rdft2.o 

C_DEPS += \
./external/fftw-3.3.1/threads/api.d \
./external/fftw-3.3.1/threads/conf.d \
./external/fftw-3.3.1/threads/ct.d \
./external/fftw-3.3.1/threads/dft-vrank-geq1.d \
./external/fftw-3.3.1/threads/f77api.d \
./external/fftw-3.3.1/threads/hc2hc.d \
./external/fftw-3.3.1/threads/openmp.d \
./external/fftw-3.3.1/threads/rdft-vrank-geq1.d \
./external/fftw-3.3.1/threads/threads.d \
./external/fftw-3.3.1/threads/vrank-geq1-rdft2.d 


# Each subdirectory must supply rules for building sources it contributes
external/fftw-3.3.1/threads/%.o: ../external/fftw-3.3.1/threads/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


