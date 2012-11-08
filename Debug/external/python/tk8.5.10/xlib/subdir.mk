################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/tk8.5.10/xlib/xcolors.c \
../external/python/tk8.5.10/xlib/xdraw.c \
../external/python/tk8.5.10/xlib/xgc.c \
../external/python/tk8.5.10/xlib/ximage.c \
../external/python/tk8.5.10/xlib/xutil.c 

OBJS += \
./external/python/tk8.5.10/xlib/xcolors.o \
./external/python/tk8.5.10/xlib/xdraw.o \
./external/python/tk8.5.10/xlib/xgc.o \
./external/python/tk8.5.10/xlib/ximage.o \
./external/python/tk8.5.10/xlib/xutil.o 

C_DEPS += \
./external/python/tk8.5.10/xlib/xcolors.d \
./external/python/tk8.5.10/xlib/xdraw.d \
./external/python/tk8.5.10/xlib/xgc.d \
./external/python/tk8.5.10/xlib/ximage.d \
./external/python/tk8.5.10/xlib/xutil.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/tk8.5.10/xlib/%.o: ../external/python/tk8.5.10/xlib/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


