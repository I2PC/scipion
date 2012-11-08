################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/numpy-1.6.1/numpy/core/src/umath/ufunc_object.c \
../external/python/numpy-1.6.1/numpy/core/src/umath/umathmodule_onefile.c 

OBJS += \
./external/python/numpy-1.6.1/numpy/core/src/umath/ufunc_object.o \
./external/python/numpy-1.6.1/numpy/core/src/umath/umathmodule_onefile.o 

C_DEPS += \
./external/python/numpy-1.6.1/numpy/core/src/umath/ufunc_object.d \
./external/python/numpy-1.6.1/numpy/core/src/umath/umathmodule_onefile.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/numpy-1.6.1/numpy/core/src/umath/%.o: ../external/python/numpy-1.6.1/numpy/core/src/umath/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


