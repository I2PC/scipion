################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/ffi.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/ffi_darwin.c 

S_UPPER_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/aix.S \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/aix_closure.S \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/darwin.S \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/darwin_closure.S \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/linux64.S \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/linux64_closure.S \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/ppc_closure.S \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/sysv.S 

OBJS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/aix.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/aix_closure.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/darwin.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/darwin_closure.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/ffi.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/ffi_darwin.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/linux64.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/linux64_closure.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/ppc_closure.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/sysv.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/ffi.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/ffi_darwin.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/%.S
	@echo 'Building file: $<'
	@echo 'Invoking: GCC Assembler'
	as  -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/powerpc/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


