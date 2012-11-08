################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/ffi.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/ffi64.c 

S_UPPER_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/darwin.S \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/darwin64.S \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/freebsd.S \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/sysv.S \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/unix64.S \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/win32.S \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/win64.S 

OBJS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/darwin.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/darwin64.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/ffi.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/ffi64.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/freebsd.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/sysv.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/unix64.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/win32.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/win64.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/ffi.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/ffi64.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/%.S
	@echo 'Building file: $<'
	@echo 'Invoking: GCC Assembler'
	as  -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/x86/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


