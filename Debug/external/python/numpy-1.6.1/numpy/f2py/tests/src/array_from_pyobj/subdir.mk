################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/numpy-1.6.1/numpy/f2py/tests/src/array_from_pyobj/wrapmodule.c 

OBJS += \
./external/python/numpy-1.6.1/numpy/f2py/tests/src/array_from_pyobj/wrapmodule.o 

C_DEPS += \
./external/python/numpy-1.6.1/numpy/f2py/tests/src/array_from_pyobj/wrapmodule.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/numpy-1.6.1/numpy/f2py/tests/src/array_from_pyobj/%.o: ../external/python/numpy-1.6.1/numpy/f2py/tests/src/array_from_pyobj/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


