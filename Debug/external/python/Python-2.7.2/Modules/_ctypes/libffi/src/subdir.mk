################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/closures.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/debug.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/dlmalloc.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/java_raw_api.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/prep_cif.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/raw_api.c \
../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/types.c 

OBJS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/closures.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/debug.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/dlmalloc.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/java_raw_api.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/prep_cif.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/raw_api.o \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/types.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/closures.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/debug.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/dlmalloc.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/java_raw_api.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/prep_cif.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/raw_api.d \
./external/python/Python-2.7.2/Modules/_ctypes/libffi/src/types.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/_ctypes/libffi/src/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/libffi/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


