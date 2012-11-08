################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/_io/_iomodule.c \
../external/python/Python-2.7.2/Modules/_io/bufferedio.c \
../external/python/Python-2.7.2/Modules/_io/bytesio.c \
../external/python/Python-2.7.2/Modules/_io/fileio.c \
../external/python/Python-2.7.2/Modules/_io/iobase.c \
../external/python/Python-2.7.2/Modules/_io/stringio.c \
../external/python/Python-2.7.2/Modules/_io/textio.c 

OBJS += \
./external/python/Python-2.7.2/Modules/_io/_iomodule.o \
./external/python/Python-2.7.2/Modules/_io/bufferedio.o \
./external/python/Python-2.7.2/Modules/_io/bytesio.o \
./external/python/Python-2.7.2/Modules/_io/fileio.o \
./external/python/Python-2.7.2/Modules/_io/iobase.o \
./external/python/Python-2.7.2/Modules/_io/stringio.o \
./external/python/Python-2.7.2/Modules/_io/textio.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/_io/_iomodule.d \
./external/python/Python-2.7.2/Modules/_io/bufferedio.d \
./external/python/Python-2.7.2/Modules/_io/bytesio.d \
./external/python/Python-2.7.2/Modules/_io/fileio.d \
./external/python/Python-2.7.2/Modules/_io/iobase.d \
./external/python/Python-2.7.2/Modules/_io/stringio.d \
./external/python/Python-2.7.2/Modules/_io/textio.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/_io/%.o: ../external/python/Python-2.7.2/Modules/_io/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


