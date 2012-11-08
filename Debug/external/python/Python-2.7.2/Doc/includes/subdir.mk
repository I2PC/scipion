################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Doc/includes/noddy.c \
../external/python/Python-2.7.2/Doc/includes/noddy2.c \
../external/python/Python-2.7.2/Doc/includes/noddy3.c \
../external/python/Python-2.7.2/Doc/includes/noddy4.c \
../external/python/Python-2.7.2/Doc/includes/run-func.c \
../external/python/Python-2.7.2/Doc/includes/shoddy.c 

OBJS += \
./external/python/Python-2.7.2/Doc/includes/noddy.o \
./external/python/Python-2.7.2/Doc/includes/noddy2.o \
./external/python/Python-2.7.2/Doc/includes/noddy3.o \
./external/python/Python-2.7.2/Doc/includes/noddy4.o \
./external/python/Python-2.7.2/Doc/includes/run-func.o \
./external/python/Python-2.7.2/Doc/includes/shoddy.o 

C_DEPS += \
./external/python/Python-2.7.2/Doc/includes/noddy.d \
./external/python/Python-2.7.2/Doc/includes/noddy2.d \
./external/python/Python-2.7.2/Doc/includes/noddy3.d \
./external/python/Python-2.7.2/Doc/includes/noddy4.d \
./external/python/Python-2.7.2/Doc/includes/run-func.d \
./external/python/Python-2.7.2/Doc/includes/shoddy.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Doc/includes/%.o: ../external/python/Python-2.7.2/Doc/includes/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


