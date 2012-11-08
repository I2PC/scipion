################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../external/sqlite-3.6.23/shell.o \
../external/sqlite-3.6.23/sqlite3.o 

C_SRCS += \
../external/sqlite-3.6.23/shell.c \
../external/sqlite-3.6.23/sqlite3.c 

OBJS += \
./external/sqlite-3.6.23/shell.o \
./external/sqlite-3.6.23/sqlite3.o 

C_DEPS += \
./external/sqlite-3.6.23/shell.d \
./external/sqlite-3.6.23/sqlite3.d 


# Each subdirectory must supply rules for building sources it contributes
external/sqlite-3.6.23/%.o: ../external/sqlite-3.6.23/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


