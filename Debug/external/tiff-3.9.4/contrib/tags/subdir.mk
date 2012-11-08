################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/tiff-3.9.4/contrib/tags/listtif.c \
../external/tiff-3.9.4/contrib/tags/maketif.c \
../external/tiff-3.9.4/contrib/tags/xtif_dir.c 

OBJS += \
./external/tiff-3.9.4/contrib/tags/listtif.o \
./external/tiff-3.9.4/contrib/tags/maketif.o \
./external/tiff-3.9.4/contrib/tags/xtif_dir.o 

C_DEPS += \
./external/tiff-3.9.4/contrib/tags/listtif.d \
./external/tiff-3.9.4/contrib/tags/maketif.d \
./external/tiff-3.9.4/contrib/tags/xtif_dir.d 


# Each subdirectory must supply rules for building sources it contributes
external/tiff-3.9.4/contrib/tags/%.o: ../external/tiff-3.9.4/contrib/tags/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


