################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/_multiprocessing/multiprocessing.c \
../external/python/Python-2.7.2/Modules/_multiprocessing/pipe_connection.c \
../external/python/Python-2.7.2/Modules/_multiprocessing/semaphore.c \
../external/python/Python-2.7.2/Modules/_multiprocessing/socket_connection.c \
../external/python/Python-2.7.2/Modules/_multiprocessing/win32_functions.c 

OBJS += \
./external/python/Python-2.7.2/Modules/_multiprocessing/multiprocessing.o \
./external/python/Python-2.7.2/Modules/_multiprocessing/pipe_connection.o \
./external/python/Python-2.7.2/Modules/_multiprocessing/semaphore.o \
./external/python/Python-2.7.2/Modules/_multiprocessing/socket_connection.o \
./external/python/Python-2.7.2/Modules/_multiprocessing/win32_functions.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/_multiprocessing/multiprocessing.d \
./external/python/Python-2.7.2/Modules/_multiprocessing/pipe_connection.d \
./external/python/Python-2.7.2/Modules/_multiprocessing/semaphore.d \
./external/python/Python-2.7.2/Modules/_multiprocessing/socket_connection.d \
./external/python/Python-2.7.2/Modules/_multiprocessing/win32_functions.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/_multiprocessing/%.o: ../external/python/Python-2.7.2/Modules/_multiprocessing/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


