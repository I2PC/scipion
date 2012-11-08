################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/PC/os2emx/config.c \
../external/python/Python-2.7.2/PC/os2emx/dlfcn.c \
../external/python/Python-2.7.2/PC/os2emx/dllentry.c \
../external/python/Python-2.7.2/PC/os2emx/getpathp.c \
../external/python/Python-2.7.2/PC/os2emx/pythonpm.c 

OBJS += \
./external/python/Python-2.7.2/PC/os2emx/config.o \
./external/python/Python-2.7.2/PC/os2emx/dlfcn.o \
./external/python/Python-2.7.2/PC/os2emx/dllentry.o \
./external/python/Python-2.7.2/PC/os2emx/getpathp.o \
./external/python/Python-2.7.2/PC/os2emx/pythonpm.o 

C_DEPS += \
./external/python/Python-2.7.2/PC/os2emx/config.d \
./external/python/Python-2.7.2/PC/os2emx/dlfcn.d \
./external/python/Python-2.7.2/PC/os2emx/dllentry.d \
./external/python/Python-2.7.2/PC/os2emx/getpathp.d \
./external/python/Python-2.7.2/PC/os2emx/pythonpm.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/PC/os2emx/%.o: ../external/python/Python-2.7.2/PC/os2emx/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


