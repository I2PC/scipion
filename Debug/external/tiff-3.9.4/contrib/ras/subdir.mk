################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/tiff-3.9.4/contrib/ras/ras2tif.c \
../external/tiff-3.9.4/contrib/ras/tif2ras.c 

OBJS += \
./external/tiff-3.9.4/contrib/ras/ras2tif.o \
./external/tiff-3.9.4/contrib/ras/tif2ras.o 

C_DEPS += \
./external/tiff-3.9.4/contrib/ras/ras2tif.d \
./external/tiff-3.9.4/contrib/ras/tif2ras.d 


# Each subdirectory must supply rules for building sources it contributes
external/tiff-3.9.4/contrib/ras/%.o: ../external/tiff-3.9.4/contrib/ras/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


