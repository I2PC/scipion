################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/ffi.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/types.c 

OBJS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/ffi.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/types.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/ffi.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/types.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


