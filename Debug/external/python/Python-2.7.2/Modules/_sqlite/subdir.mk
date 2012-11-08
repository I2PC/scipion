################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/_sqlite/cache.c \
../external/python/Python-2.7.2/Modules/_sqlite/connection.c \
../external/python/Python-2.7.2/Modules/_sqlite/cursor.c \
../external/python/Python-2.7.2/Modules/_sqlite/microprotocols.c \
../external/python/Python-2.7.2/Modules/_sqlite/module.c \
../external/python/Python-2.7.2/Modules/_sqlite/prepare_protocol.c \
../external/python/Python-2.7.2/Modules/_sqlite/row.c \
../external/python/Python-2.7.2/Modules/_sqlite/statement.c \
../external/python/Python-2.7.2/Modules/_sqlite/util.c 

OBJS += \
./external/python/Python-2.7.2/Modules/_sqlite/cache.o \
./external/python/Python-2.7.2/Modules/_sqlite/connection.o \
./external/python/Python-2.7.2/Modules/_sqlite/cursor.o \
./external/python/Python-2.7.2/Modules/_sqlite/microprotocols.o \
./external/python/Python-2.7.2/Modules/_sqlite/module.o \
./external/python/Python-2.7.2/Modules/_sqlite/prepare_protocol.o \
./external/python/Python-2.7.2/Modules/_sqlite/row.o \
./external/python/Python-2.7.2/Modules/_sqlite/statement.o \
./external/python/Python-2.7.2/Modules/_sqlite/util.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/_sqlite/cache.d \
./external/python/Python-2.7.2/Modules/_sqlite/connection.d \
./external/python/Python-2.7.2/Modules/_sqlite/cursor.d \
./external/python/Python-2.7.2/Modules/_sqlite/microprotocols.d \
./external/python/Python-2.7.2/Modules/_sqlite/module.d \
./external/python/Python-2.7.2/Modules/_sqlite/prepare_protocol.d \
./external/python/Python-2.7.2/Modules/_sqlite/row.d \
./external/python/Python-2.7.2/Modules/_sqlite/statement.d \
./external/python/Python-2.7.2/Modules/_sqlite/util.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/_sqlite/%.o: ../external/python/Python-2.7.2/Modules/_sqlite/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


