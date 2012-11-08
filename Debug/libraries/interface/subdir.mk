################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../libraries/interface/docfile.cpp \
../libraries/interface/selfile.cpp \
../libraries/interface/spider.cpp 

OBJS += \
./libraries/interface/docfile.o \
./libraries/interface/selfile.o \
./libraries/interface/spider.o 

CPP_DEPS += \
./libraries/interface/docfile.d \
./libraries/interface/selfile.d \
./libraries/interface/spider.d 


# Each subdirectory must supply rules for building sources it contributes
libraries/interface/%.o: ../libraries/interface/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


