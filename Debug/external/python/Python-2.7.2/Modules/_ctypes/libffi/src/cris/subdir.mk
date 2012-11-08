################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/cris/ffi.c 

S_UPPER_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/cris/sysv.S 

OBJS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/cris/ffi.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/cris/sysv.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/cris/ffi.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/_ctypes/libffi/src/cris/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/cris/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/python/Python-2.7.2/Modules/_ctypes/libffi/src/cris/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/cris/%.S
	@echo 'Building file: $<'
	@echo 'Invoking: GCC Assembler'
	as  -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


