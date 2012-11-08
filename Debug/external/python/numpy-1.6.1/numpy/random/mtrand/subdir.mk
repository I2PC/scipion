################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/numpy-1.6.1/numpy/random/mtrand/distributions.c \
../external/python/numpy-1.6.1/numpy/random/mtrand/initarray.c \
../external/python/numpy-1.6.1/numpy/random/mtrand/mtrand.c \
../external/python/numpy-1.6.1/numpy/random/mtrand/randomkit.c 

OBJS += \
./external/python/numpy-1.6.1/numpy/random/mtrand/distributions.o \
./external/python/numpy-1.6.1/numpy/random/mtrand/initarray.o \
./external/python/numpy-1.6.1/numpy/random/mtrand/mtrand.o \
./external/python/numpy-1.6.1/numpy/random/mtrand/randomkit.o 

C_DEPS += \
./external/python/numpy-1.6.1/numpy/random/mtrand/distributions.d \
./external/python/numpy-1.6.1/numpy/random/mtrand/initarray.d \
./external/python/numpy-1.6.1/numpy/random/mtrand/mtrand.d \
./external/python/numpy-1.6.1/numpy/random/mtrand/randomkit.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/numpy-1.6.1/numpy/random/mtrand/%.o: ../external/python/numpy-1.6.1/numpy/random/mtrand/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


