################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/tiff-3.9.4/contrib/mac-cw/mac_main.c \
../external/tiff-3.9.4/contrib/mac-cw/mkg3_main.c 

OBJS += \
./external/tiff-3.9.4/contrib/mac-cw/mac_main.o \
./external/tiff-3.9.4/contrib/mac-cw/mkg3_main.o 

C_DEPS += \
./external/tiff-3.9.4/contrib/mac-cw/mac_main.d \
./external/tiff-3.9.4/contrib/mac-cw/mkg3_main.d 


# Each subdirectory must supply rules for building sources it contributes
external/tiff-3.9.4/contrib/mac-cw/%.o: ../external/tiff-3.9.4/contrib/mac-cw/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


