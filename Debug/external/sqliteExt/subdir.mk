################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/sqliteExt/extension-functions.c 

OBJS += \
./external/sqliteExt/extension-functions.o 

C_DEPS += \
./external/sqliteExt/extension-functions.d 


# Each subdirectory must supply rules for building sources it contributes
external/sqliteExt/%.o: ../external/sqliteExt/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


