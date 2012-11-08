################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/tiff-3.9.4/test/ascii_tag.c \
../external/tiff-3.9.4/test/check_tag.c \
../external/tiff-3.9.4/test/long_tag.c \
../external/tiff-3.9.4/test/short_tag.c \
../external/tiff-3.9.4/test/strip.c \
../external/tiff-3.9.4/test/strip_rw.c \
../external/tiff-3.9.4/test/test_arrays.c 

OBJS += \
./external/tiff-3.9.4/test/ascii_tag.o \
./external/tiff-3.9.4/test/check_tag.o \
./external/tiff-3.9.4/test/long_tag.o \
./external/tiff-3.9.4/test/short_tag.o \
./external/tiff-3.9.4/test/strip.o \
./external/tiff-3.9.4/test/strip_rw.o \
./external/tiff-3.9.4/test/test_arrays.o 

C_DEPS += \
./external/tiff-3.9.4/test/ascii_tag.d \
./external/tiff-3.9.4/test/check_tag.d \
./external/tiff-3.9.4/test/long_tag.d \
./external/tiff-3.9.4/test/short_tag.d \
./external/tiff-3.9.4/test/strip.d \
./external/tiff-3.9.4/test/strip_rw.d \
./external/tiff-3.9.4/test/test_arrays.d 


# Each subdirectory must supply rules for building sources it contributes
external/tiff-3.9.4/test/%.o: ../external/tiff-3.9.4/test/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


