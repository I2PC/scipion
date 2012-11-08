################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../libraries/bindings/python/python_constants.cpp \
../libraries/bindings/python/python_filename.cpp \
../libraries/bindings/python/python_image.cpp \
../libraries/bindings/python/python_metadata.cpp \
../libraries/bindings/python/python_program.cpp \
../libraries/bindings/python/python_symmetry.cpp \
../libraries/bindings/python/xmippmodule.cpp 

OBJS += \
./libraries/bindings/python/python_constants.o \
./libraries/bindings/python/python_filename.o \
./libraries/bindings/python/python_image.o \
./libraries/bindings/python/python_metadata.o \
./libraries/bindings/python/python_program.o \
./libraries/bindings/python/python_symmetry.o \
./libraries/bindings/python/xmippmodule.o 

CPP_DEPS += \
./libraries/bindings/python/python_constants.d \
./libraries/bindings/python/python_filename.d \
./libraries/bindings/python/python_image.d \
./libraries/bindings/python/python_metadata.d \
./libraries/bindings/python/python_program.d \
./libraries/bindings/python/python_symmetry.d \
./libraries/bindings/python/xmippmodule.d 


# Each subdirectory must supply rules for building sources it contributes
libraries/bindings/python/%.o: ../libraries/bindings/python/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


