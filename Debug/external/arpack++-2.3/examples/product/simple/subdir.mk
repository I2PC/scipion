################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../external/arpack++-2.3/examples/product/simple/symsimp.cc 

OBJS += \
./external/arpack++-2.3/examples/product/simple/symsimp.o 

CC_DEPS += \
./external/arpack++-2.3/examples/product/simple/symsimp.d 


# Each subdirectory must supply rules for building sources it contributes
external/arpack++-2.3/examples/product/simple/%.o: ../external/arpack++-2.3/examples/product/simple/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


