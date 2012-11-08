################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/matplotlib-1.1.0/CXX/cxxextensions.c 

CXX_SRCS += \
../external/python/matplotlib-1.1.0/CXX/IndirectPythonInterface.cxx \
../external/python/matplotlib-1.1.0/CXX/cxx_extensions.cxx \
../external/python/matplotlib-1.1.0/CXX/cxxsupport.cxx 

OBJS += \
./external/python/matplotlib-1.1.0/CXX/IndirectPythonInterface.o \
./external/python/matplotlib-1.1.0/CXX/cxx_extensions.o \
./external/python/matplotlib-1.1.0/CXX/cxxextensions.o \
./external/python/matplotlib-1.1.0/CXX/cxxsupport.o 

C_DEPS += \
./external/python/matplotlib-1.1.0/CXX/cxxextensions.d 

CXX_DEPS += \
./external/python/matplotlib-1.1.0/CXX/IndirectPythonInterface.d \
./external/python/matplotlib-1.1.0/CXX/cxx_extensions.d \
./external/python/matplotlib-1.1.0/CXX/cxxsupport.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/matplotlib-1.1.0/CXX/%.o: ../external/python/matplotlib-1.1.0/CXX/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/python/matplotlib-1.1.0/CXX/%.o: ../external/python/matplotlib-1.1.0/CXX/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


