################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/sh/ffi.c 

S_UPPER_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/sh/sysv.S 

OBJS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/sh/ffi.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/sh/sysv.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/sh/ffi.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/_ctypes/libffi/src/sh/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/sh/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/python/Python-2.7.2/Modules/_ctypes/libffi/src/sh/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/sh/%.S
	@echo 'Building file: $<'
	@echo 'Invoking: GCC Assembler'
	as  -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


