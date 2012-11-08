################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Mac/Modules/ColorPickermodule.c \
../external/python/Python-2.7.2/Mac/Modules/MacOS.c \
../external/python/Python-2.7.2/Mac/Modules/Nav.c \
../external/python/Python-2.7.2/Mac/Modules/OSATerminology.c \
../external/python/Python-2.7.2/Mac/Modules/_scproxy.c \
../external/python/Python-2.7.2/Mac/Modules/autoGIL.c \
../external/python/Python-2.7.2/Mac/Modules/gestaltmodule.c \
../external/python/Python-2.7.2/Mac/Modules/icgluemodule.c 

OBJS += \
./external/python/Python-2.7.2/Mac/Modules/ColorPickermodule.o \
./external/python/Python-2.7.2/Mac/Modules/MacOS.o \
./external/python/Python-2.7.2/Mac/Modules/Nav.o \
./external/python/Python-2.7.2/Mac/Modules/OSATerminology.o \
./external/python/Python-2.7.2/Mac/Modules/_scproxy.o \
./external/python/Python-2.7.2/Mac/Modules/autoGIL.o \
./external/python/Python-2.7.2/Mac/Modules/gestaltmodule.o \
./external/python/Python-2.7.2/Mac/Modules/icgluemodule.o 

C_DEPS += \
./external/python/Python-2.7.2/Mac/Modules/ColorPickermodule.d \
./external/python/Python-2.7.2/Mac/Modules/MacOS.d \
./external/python/Python-2.7.2/Mac/Modules/Nav.d \
./external/python/Python-2.7.2/Mac/Modules/OSATerminology.d \
./external/python/Python-2.7.2/Mac/Modules/_scproxy.d \
./external/python/Python-2.7.2/Mac/Modules/autoGIL.d \
./external/python/Python-2.7.2/Mac/Modules/gestaltmodule.d \
./external/python/Python-2.7.2/Mac/Modules/icgluemodule.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Mac/Modules/%.o: ../external/python/Python-2.7.2/Mac/Modules/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


