################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/PC/os2vacpp/config.c \
../external/python/Python-2.7.2/PC/os2vacpp/getpathp.c 

OBJS += \
./external/python/Python-2.7.2/PC/os2vacpp/config.o \
./external/python/Python-2.7.2/PC/os2vacpp/getpathp.o 

C_DEPS += \
./external/python/Python-2.7.2/PC/os2vacpp/config.d \
./external/python/Python-2.7.2/PC/os2vacpp/getpathp.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/PC/os2vacpp/%.o: ../external/python/Python-2.7.2/PC/os2vacpp/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


