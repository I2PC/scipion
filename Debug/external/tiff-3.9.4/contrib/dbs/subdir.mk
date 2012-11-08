################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../external/tiff-3.9.4/contrib/dbs/tiff-bi.o \
../external/tiff-3.9.4/contrib/dbs/tiff-grayscale.o \
../external/tiff-3.9.4/contrib/dbs/tiff-palette.o \
../external/tiff-3.9.4/contrib/dbs/tiff-rgb.o 

C_SRCS += \
../external/tiff-3.9.4/contrib/dbs/tiff-bi.c \
../external/tiff-3.9.4/contrib/dbs/tiff-grayscale.c \
../external/tiff-3.9.4/contrib/dbs/tiff-palette.c \
../external/tiff-3.9.4/contrib/dbs/tiff-rgb.c 

OBJS += \
./external/tiff-3.9.4/contrib/dbs/tiff-bi.o \
./external/tiff-3.9.4/contrib/dbs/tiff-grayscale.o \
./external/tiff-3.9.4/contrib/dbs/tiff-palette.o \
./external/tiff-3.9.4/contrib/dbs/tiff-rgb.o 

C_DEPS += \
./external/tiff-3.9.4/contrib/dbs/tiff-bi.d \
./external/tiff-3.9.4/contrib/dbs/tiff-grayscale.d \
./external/tiff-3.9.4/contrib/dbs/tiff-palette.d \
./external/tiff-3.9.4/contrib/dbs/tiff-rgb.d 


# Each subdirectory must supply rules for building sources it contributes
external/tiff-3.9.4/contrib/dbs/%.o: ../external/tiff-3.9.4/contrib/dbs/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


