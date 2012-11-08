################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../external/tiff-3.9.4/contrib/iptcutil/iptcutil.o 

C_SRCS += \
../external/tiff-3.9.4/contrib/iptcutil/iptcutil.c 

OBJS += \
./external/tiff-3.9.4/contrib/iptcutil/iptcutil.o 

C_DEPS += \
./external/tiff-3.9.4/contrib/iptcutil/iptcutil.d 


# Each subdirectory must supply rules for building sources it contributes
external/tiff-3.9.4/contrib/iptcutil/%.o: ../external/tiff-3.9.4/contrib/iptcutil/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


