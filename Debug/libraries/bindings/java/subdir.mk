################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../libraries/bindings/java/xmipp_Aux.cpp \
../libraries/bindings/java/xmipp_CTFDescription.cpp \
../libraries/bindings/java/xmipp_ExceptionsHandler.cpp \
../libraries/bindings/java/xmipp_Filename.cpp \
../libraries/bindings/java/xmipp_ImageGeneric.cpp \
../libraries/bindings/java/xmipp_MetaData.cpp \
../libraries/bindings/java/xmipp_ProgTomographAlignment.cpp \
../libraries/bindings/java/xmipp_Program.cpp \
../libraries/bindings/java/xmipp_Projection.cpp \
../libraries/bindings/java/xmipp_TiltPairAligner.cpp 

OBJS += \
./libraries/bindings/java/xmipp_Aux.o \
./libraries/bindings/java/xmipp_CTFDescription.o \
./libraries/bindings/java/xmipp_ExceptionsHandler.o \
./libraries/bindings/java/xmipp_Filename.o \
./libraries/bindings/java/xmipp_ImageGeneric.o \
./libraries/bindings/java/xmipp_MetaData.o \
./libraries/bindings/java/xmipp_ProgTomographAlignment.o \
./libraries/bindings/java/xmipp_Program.o \
./libraries/bindings/java/xmipp_Projection.o \
./libraries/bindings/java/xmipp_TiltPairAligner.o 

CPP_DEPS += \
./libraries/bindings/java/xmipp_Aux.d \
./libraries/bindings/java/xmipp_CTFDescription.d \
./libraries/bindings/java/xmipp_ExceptionsHandler.d \
./libraries/bindings/java/xmipp_Filename.d \
./libraries/bindings/java/xmipp_ImageGeneric.d \
./libraries/bindings/java/xmipp_MetaData.d \
./libraries/bindings/java/xmipp_ProgTomographAlignment.d \
./libraries/bindings/java/xmipp_Program.d \
./libraries/bindings/java/xmipp_Projection.d \
./libraries/bindings/java/xmipp_TiltPairAligner.d 


# Each subdirectory must supply rules for building sources it contributes
libraries/bindings/java/%.o: ../libraries/bindings/java/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


