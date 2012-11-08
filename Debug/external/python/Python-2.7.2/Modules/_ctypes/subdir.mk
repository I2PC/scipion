################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/_ctypes/_ctypes.c \
../external/python/Python-2.7.2/Modules/_ctypes/_ctypes_test.c \
../external/python/Python-2.7.2/Modules/_ctypes/callbacks.c \
../external/python/Python-2.7.2/Modules/_ctypes/callproc.c \
../external/python/Python-2.7.2/Modules/_ctypes/cfield.c \
../external/python/Python-2.7.2/Modules/_ctypes/malloc_closure.c \
../external/python/Python-2.7.2/Modules/_ctypes/stgdict.c 

OBJS += \
./external/python/Python-2.7.2/Modules/_ctypes/_ctypes.o \
./external/python/Python-2.7.2/Modules/_ctypes/_ctypes_test.o \
./external/python/Python-2.7.2/Modules/_ctypes/callbacks.o \
./external/python/Python-2.7.2/Modules/_ctypes/callproc.o \
./external/python/Python-2.7.2/Modules/_ctypes/cfield.o \
./external/python/Python-2.7.2/Modules/_ctypes/malloc_closure.o \
./external/python/Python-2.7.2/Modules/_ctypes/stgdict.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/_ctypes/_ctypes.d \
./external/python/Python-2.7.2/Modules/_ctypes/_ctypes_test.d \
./external/python/Python-2.7.2/Modules/_ctypes/callbacks.d \
./external/python/Python-2.7.2/Modules/_ctypes/callproc.d \
./external/python/Python-2.7.2/Modules/_ctypes/cfield.d \
./external/python/Python-2.7.2/Modules/_ctypes/malloc_closure.d \
./external/python/Python-2.7.2/Modules/_ctypes/stgdict.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/_ctypes/%.o: ../external/python/Python-2.7.2/Modules/_ctypes/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


