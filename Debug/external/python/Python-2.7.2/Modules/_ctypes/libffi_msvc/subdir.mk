################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/ffi.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/prep_cif.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/types.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/win32.c 

ASM_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/win64.asm 

OBJS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/ffi.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/prep_cif.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/types.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/win32.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/win64.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/ffi.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/prep_cif.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/types.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/win32.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi_msvc/%.asm
	@echo 'Building file: $<'
	@echo 'Invoking: GCC Assembler'
	as  -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


