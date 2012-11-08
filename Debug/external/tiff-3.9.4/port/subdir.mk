################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../external/tiff-3.9.4/port/dummy.o 

C_SRCS += \
../external/tiff-3.9.4/port/dummy.c \
../external/tiff-3.9.4/port/getopt.c \
../external/tiff-3.9.4/port/lfind.c \
../external/tiff-3.9.4/port/strcasecmp.c \
../external/tiff-3.9.4/port/strtoul.c 

OBJS += \
./external/tiff-3.9.4/port/dummy.o \
./external/tiff-3.9.4/port/getopt.o \
./external/tiff-3.9.4/port/lfind.o \
./external/tiff-3.9.4/port/strcasecmp.o \
./external/tiff-3.9.4/port/strtoul.o 

C_DEPS += \
./external/tiff-3.9.4/port/dummy.d \
./external/tiff-3.9.4/port/getopt.d \
./external/tiff-3.9.4/port/lfind.d \
./external/tiff-3.9.4/port/strcasecmp.d \
./external/tiff-3.9.4/port/strtoul.d 


# Each subdirectory must supply rules for building sources it contributes
external/tiff-3.9.4/port/%.o: ../external/tiff-3.9.4/port/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


