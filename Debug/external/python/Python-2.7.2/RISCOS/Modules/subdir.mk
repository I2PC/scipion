################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/RISCOS/Modules/config.c \
../external/python/Python-2.7.2/RISCOS/Modules/drawfmodule.c \
../external/python/Python-2.7.2/RISCOS/Modules/getpath_riscos.c \
../external/python/Python-2.7.2/RISCOS/Modules/riscosmodule.c \
../external/python/Python-2.7.2/RISCOS/Modules/swimodule.c 

OBJS += \
./external/python/Python-2.7.2/RISCOS/Modules/config.o \
./external/python/Python-2.7.2/RISCOS/Modules/drawfmodule.o \
./external/python/Python-2.7.2/RISCOS/Modules/getpath_riscos.o \
./external/python/Python-2.7.2/RISCOS/Modules/riscosmodule.o \
./external/python/Python-2.7.2/RISCOS/Modules/swimodule.o 

C_DEPS += \
./external/python/Python-2.7.2/RISCOS/Modules/config.d \
./external/python/Python-2.7.2/RISCOS/Modules/drawfmodule.d \
./external/python/Python-2.7.2/RISCOS/Modules/getpath_riscos.d \
./external/python/Python-2.7.2/RISCOS/Modules/riscosmodule.d \
./external/python/Python-2.7.2/RISCOS/Modules/swimodule.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/RISCOS/Modules/%.o: ../external/python/Python-2.7.2/RISCOS/Modules/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


