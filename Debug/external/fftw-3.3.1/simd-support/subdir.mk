################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../external/fftw-3.3.1/simd-support/altivec.o \
../external/fftw-3.3.1/simd-support/avx.o \
../external/fftw-3.3.1/simd-support/libsimd_sse2_nonportable_la-sse2-nonportable.o \
../external/fftw-3.3.1/simd-support/neon.o \
../external/fftw-3.3.1/simd-support/sse2.o \
../external/fftw-3.3.1/simd-support/taint.o 

C_SRCS += \
../external/fftw-3.3.1/simd-support/altivec.c \
../external/fftw-3.3.1/simd-support/avx.c \
../external/fftw-3.3.1/simd-support/neon.c \
../external/fftw-3.3.1/simd-support/sse2-nonportable.c \
../external/fftw-3.3.1/simd-support/sse2.c \
../external/fftw-3.3.1/simd-support/taint.c 

OBJS += \
./external/fftw-3.3.1/simd-support/altivec.o \
./external/fftw-3.3.1/simd-support/avx.o \
./external/fftw-3.3.1/simd-support/neon.o \
./external/fftw-3.3.1/simd-support/sse2-nonportable.o \
./external/fftw-3.3.1/simd-support/sse2.o \
./external/fftw-3.3.1/simd-support/taint.o 

C_DEPS += \
./external/fftw-3.3.1/simd-support/altivec.d \
./external/fftw-3.3.1/simd-support/avx.d \
./external/fftw-3.3.1/simd-support/neon.d \
./external/fftw-3.3.1/simd-support/sse2-nonportable.d \
./external/fftw-3.3.1/simd-support/sse2.d \
./external/fftw-3.3.1/simd-support/taint.d 


# Each subdirectory must supply rules for building sources it contributes
external/fftw-3.3.1/simd-support/%.o: ../external/fftw-3.3.1/simd-support/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


