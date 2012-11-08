################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/matplotlib-1.1.0/CXX/Python3/cxxextensions.c 

CXX_SRCS += \
../external/python/matplotlib-1.1.0/CXX/Python3/IndirectPythonInterface.cxx \
../external/python/matplotlib-1.1.0/CXX/Python3/cxx_extensions.cxx \
../external/python/matplotlib-1.1.0/CXX/Python3/cxxsupport.cxx 

OBJS += \
./external/python/matplotlib-1.1.0/CXX/Python3/IndirectPythonInterface.o \
./external/python/matplotlib-1.1.0/CXX/Python3/cxx_extensions.o \
./external/python/matplotlib-1.1.0/CXX/Python3/cxxextensions.o \
./external/python/matplotlib-1.1.0/CXX/Python3/cxxsupport.o 

C_DEPS += \
./external/python/matplotlib-1.1.0/CXX/Python3/cxxextensions.d 

CXX_DEPS += \
./external/python/matplotlib-1.1.0/CXX/Python3/IndirectPythonInterface.d \
./external/python/matplotlib-1.1.0/CXX/Python3/cxx_extensions.d \
./external/python/matplotlib-1.1.0/CXX/Python3/cxxsupport.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/matplotlib-1.1.0/CXX/Python3/%.o: ../external/python/matplotlib-1.1.0/CXX/Python3/%.cxx
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/python/matplotlib-1.1.0/CXX/Python3/%.o: ../external/python/matplotlib-1.1.0/CXX/Python3/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


