################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/debug.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/ffi.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/prep_cif.c 

ASM_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/sysv.asm 

OBJS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/debug.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/ffi.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/prep_cif.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/sysv.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/debug.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/ffi.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/prep_cif.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi_arm_wince/%.asm
	@echo 'Building file: $<'
	@echo 'Invoking: GCC Assembler'
	as  -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


