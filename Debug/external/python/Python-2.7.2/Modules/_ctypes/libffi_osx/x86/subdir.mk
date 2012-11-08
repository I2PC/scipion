################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/x86/x86-ffi64.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/x86/x86-ffi_darwin.c 

S_UPPER_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/x86/darwin64.S \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/x86/x86-darwin.S 

OBJS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/x86/darwin64.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/x86/x86-darwin.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/x86/x86-ffi64.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/x86/x86-ffi_darwin.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/x86/x86-ffi64.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/x86/x86-ffi_darwin.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/x86/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/x86/%.S
	@echo 'Building file: $<'
	@echo 'Invoking: GCC Assembler'
	as  -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/x86/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi_osx/x86/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


