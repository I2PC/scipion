################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/tiff-3.9.4/contrib/mfs/mfs_file.c 

OBJS += \
./external/tiff-3.9.4/contrib/mfs/mfs_file.o 

C_DEPS += \
./external/tiff-3.9.4/contrib/mfs/mfs_file.d 


# Each subdirectory must supply rules for building sources it contributes
external/tiff-3.9.4/contrib/mfs/%.o: ../external/tiff-3.9.4/contrib/mfs/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


