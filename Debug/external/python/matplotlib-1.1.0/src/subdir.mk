################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../external/python/matplotlib-1.1.0/src/_backend_agg.cpp \
../external/python/matplotlib-1.1.0/src/_gtkagg.cpp \
../external/python/matplotlib-1.1.0/src/_image.cpp \
../external/python/matplotlib-1.1.0/src/_path.cpp \
../external/python/matplotlib-1.1.0/src/_png.cpp \
../external/python/matplotlib-1.1.0/src/_tkagg.cpp \
../external/python/matplotlib-1.1.0/src/_ttconv.cpp \
../external/python/matplotlib-1.1.0/src/_windowing.cpp \
../external/python/matplotlib-1.1.0/src/agg_py_transforms.cpp \
../external/python/matplotlib-1.1.0/src/ft2font.cpp \
../external/python/matplotlib-1.1.0/src/mplutils.cpp \
../external/python/matplotlib-1.1.0/src/path_cleanup.cpp 

C_SRCS += \
../external/python/matplotlib-1.1.0/src/_backend_gdk.c \
../external/python/matplotlib-1.1.0/src/backend_gdk.c \
../external/python/matplotlib-1.1.0/src/cntr.c \
../external/python/matplotlib-1.1.0/src/nxutils.c 

OBJS += \
./external/python/matplotlib-1.1.0/src/_backend_agg.o \
./external/python/matplotlib-1.1.0/src/_backend_gdk.o \
./external/python/matplotlib-1.1.0/src/_gtkagg.o \
./external/python/matplotlib-1.1.0/src/_image.o \
./external/python/matplotlib-1.1.0/src/_path.o \
./external/python/matplotlib-1.1.0/src/_png.o \
./external/python/matplotlib-1.1.0/src/_tkagg.o \
./external/python/matplotlib-1.1.0/src/_ttconv.o \
./external/python/matplotlib-1.1.0/src/_windowing.o \
./external/python/matplotlib-1.1.0/src/agg_py_transforms.o \
./external/python/matplotlib-1.1.0/src/backend_gdk.o \
./external/python/matplotlib-1.1.0/src/cntr.o \
./external/python/matplotlib-1.1.0/src/ft2font.o \
./external/python/matplotlib-1.1.0/src/mplutils.o \
./external/python/matplotlib-1.1.0/src/nxutils.o \
./external/python/matplotlib-1.1.0/src/path_cleanup.o 

C_DEPS += \
./external/python/matplotlib-1.1.0/src/_backend_gdk.d \
./external/python/matplotlib-1.1.0/src/backend_gdk.d \
./external/python/matplotlib-1.1.0/src/cntr.d \
./external/python/matplotlib-1.1.0/src/nxutils.d 

CPP_DEPS += \
./external/python/matplotlib-1.1.0/src/_backend_agg.d \
./external/python/matplotlib-1.1.0/src/_gtkagg.d \
./external/python/matplotlib-1.1.0/src/_image.d \
./external/python/matplotlib-1.1.0/src/_path.d \
./external/python/matplotlib-1.1.0/src/_png.d \
./external/python/matplotlib-1.1.0/src/_tkagg.d \
./external/python/matplotlib-1.1.0/src/_ttconv.d \
./external/python/matplotlib-1.1.0/src/_windowing.d \
./external/python/matplotlib-1.1.0/src/agg_py_transforms.d \
./external/python/matplotlib-1.1.0/src/ft2font.d \
./external/python/matplotlib-1.1.0/src/mplutils.d \
./external/python/matplotlib-1.1.0/src/path_cleanup.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/matplotlib-1.1.0/src/%.o: ../external/python/matplotlib-1.1.0/src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

external/python/matplotlib-1.1.0/src/%.o: ../external/python/matplotlib-1.1.0/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


