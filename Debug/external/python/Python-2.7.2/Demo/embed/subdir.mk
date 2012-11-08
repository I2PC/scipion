################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Demo/embed/demo.c \
../external/python/Python-2.7.2/Demo/embed/importexc.c \
../external/python/Python-2.7.2/Demo/embed/loop.c 

OBJS += \
./external/python/Python-2.7.2/Demo/embed/demo.o \
./external/python/Python-2.7.2/Demo/embed/importexc.o \
./external/python/Python-2.7.2/Demo/embed/loop.o 

C_DEPS += \
./external/python/Python-2.7.2/Demo/embed/demo.d \
./external/python/Python-2.7.2/Demo/embed/importexc.d \
./external/python/Python-2.7.2/Demo/embed/loop.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Demo/embed/%.o: ../external/python/Python-2.7.2/Demo/embed/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


