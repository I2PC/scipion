################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/tcl8.5.10/win/cat.c \
../external/python/tcl8.5.10/win/nmakehlp.c \
../external/python/tcl8.5.10/win/stub16.c \
../external/python/tcl8.5.10/win/tclAppInit.c \
../external/python/tcl8.5.10/win/tclWin32Dll.c \
../external/python/tcl8.5.10/win/tclWinChan.c \
../external/python/tcl8.5.10/win/tclWinConsole.c \
../external/python/tcl8.5.10/win/tclWinDde.c \
../external/python/tcl8.5.10/win/tclWinError.c \
../external/python/tcl8.5.10/win/tclWinFCmd.c \
../external/python/tcl8.5.10/win/tclWinFile.c \
../external/python/tcl8.5.10/win/tclWinInit.c \
../external/python/tcl8.5.10/win/tclWinLoad.c \
../external/python/tcl8.5.10/win/tclWinNotify.c \
../external/python/tcl8.5.10/win/tclWinPipe.c \
../external/python/tcl8.5.10/win/tclWinReg.c \
../external/python/tcl8.5.10/win/tclWinSerial.c \
../external/python/tcl8.5.10/win/tclWinSock.c \
../external/python/tcl8.5.10/win/tclWinTest.c \
../external/python/tcl8.5.10/win/tclWinThrd.c \
../external/python/tcl8.5.10/win/tclWinTime.c 

OBJS += \
./external/python/tcl8.5.10/win/cat.o \
./external/python/tcl8.5.10/win/nmakehlp.o \
./external/python/tcl8.5.10/win/stub16.o \
./external/python/tcl8.5.10/win/tclAppInit.o \
./external/python/tcl8.5.10/win/tclWin32Dll.o \
./external/python/tcl8.5.10/win/tclWinChan.o \
./external/python/tcl8.5.10/win/tclWinConsole.o \
./external/python/tcl8.5.10/win/tclWinDde.o \
./external/python/tcl8.5.10/win/tclWinError.o \
./external/python/tcl8.5.10/win/tclWinFCmd.o \
./external/python/tcl8.5.10/win/tclWinFile.o \
./external/python/tcl8.5.10/win/tclWinInit.o \
./external/python/tcl8.5.10/win/tclWinLoad.o \
./external/python/tcl8.5.10/win/tclWinNotify.o \
./external/python/tcl8.5.10/win/tclWinPipe.o \
./external/python/tcl8.5.10/win/tclWinReg.o \
./external/python/tcl8.5.10/win/tclWinSerial.o \
./external/python/tcl8.5.10/win/tclWinSock.o \
./external/python/tcl8.5.10/win/tclWinTest.o \
./external/python/tcl8.5.10/win/tclWinThrd.o \
./external/python/tcl8.5.10/win/tclWinTime.o 

C_DEPS += \
./external/python/tcl8.5.10/win/cat.d \
./external/python/tcl8.5.10/win/nmakehlp.d \
./external/python/tcl8.5.10/win/stub16.d \
./external/python/tcl8.5.10/win/tclAppInit.d \
./external/python/tcl8.5.10/win/tclWin32Dll.d \
./external/python/tcl8.5.10/win/tclWinChan.d \
./external/python/tcl8.5.10/win/tclWinConsole.d \
./external/python/tcl8.5.10/win/tclWinDde.d \
./external/python/tcl8.5.10/win/tclWinError.d \
./external/python/tcl8.5.10/win/tclWinFCmd.d \
./external/python/tcl8.5.10/win/tclWinFile.d \
./external/python/tcl8.5.10/win/tclWinInit.d \
./external/python/tcl8.5.10/win/tclWinLoad.d \
./external/python/tcl8.5.10/win/tclWinNotify.d \
./external/python/tcl8.5.10/win/tclWinPipe.d \
./external/python/tcl8.5.10/win/tclWinReg.d \
./external/python/tcl8.5.10/win/tclWinSerial.d \
./external/python/tcl8.5.10/win/tclWinSock.d \
./external/python/tcl8.5.10/win/tclWinTest.d \
./external/python/tcl8.5.10/win/tclWinThrd.d \
./external/python/tcl8.5.10/win/tclWinTime.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/tcl8.5.10/win/%.o: ../external/python/tcl8.5.10/win/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


