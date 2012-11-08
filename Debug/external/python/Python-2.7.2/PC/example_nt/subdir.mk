################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/PC/example_nt/example.c 

OBJS += \
./external/python/Python-2.7.2/PC/example_nt/example.o 

C_DEPS += \
./external/python/Python-2.7.2/PC/example_nt/example.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/PC/example_nt/%.o: ../external/python/Python-2.7.2/PC/example_nt/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


