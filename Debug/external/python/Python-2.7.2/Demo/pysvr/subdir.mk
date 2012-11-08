################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Demo/pysvr/pysvr.c 

OBJS += \
./external/python/Python-2.7.2/Demo/pysvr/pysvr.o 

C_DEPS += \
./external/python/Python-2.7.2/Demo/pysvr/pysvr.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Demo/pysvr/%.o: ../external/python/Python-2.7.2/Demo/pysvr/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


