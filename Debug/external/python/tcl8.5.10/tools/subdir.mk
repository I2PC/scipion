################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/tcl8.5.10/tools/man2tcl.c 

OBJS += \
./external/python/tcl8.5.10/tools/man2tcl.o 

C_DEPS += \
./external/python/tcl8.5.10/tools/man2tcl.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/tcl8.5.10/tools/%.o: ../external/python/tcl8.5.10/tools/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


