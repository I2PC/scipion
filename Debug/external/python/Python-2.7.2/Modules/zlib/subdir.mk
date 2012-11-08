################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/zlib/adler32.c \
../external/python/Python-2.7.2/Modules/zlib/compress.c \
../external/python/Python-2.7.2/Modules/zlib/crc32.c \
../external/python/Python-2.7.2/Modules/zlib/deflate.c \
../external/python/Python-2.7.2/Modules/zlib/example.c \
../external/python/Python-2.7.2/Modules/zlib/gzio.c \
../external/python/Python-2.7.2/Modules/zlib/infback.c \
../external/python/Python-2.7.2/Modules/zlib/inffast.c \
../external/python/Python-2.7.2/Modules/zlib/inflate.c \
../external/python/Python-2.7.2/Modules/zlib/inftrees.c \
../external/python/Python-2.7.2/Modules/zlib/minigzip.c \
../external/python/Python-2.7.2/Modules/zlib/trees.c \
../external/python/Python-2.7.2/Modules/zlib/uncompr.c \
../external/python/Python-2.7.2/Modules/zlib/zutil.c 

OBJS += \
./external/python/Python-2.7.2/Modules/zlib/adler32.o \
./external/python/Python-2.7.2/Modules/zlib/compress.o \
./external/python/Python-2.7.2/Modules/zlib/crc32.o \
./external/python/Python-2.7.2/Modules/zlib/deflate.o \
./external/python/Python-2.7.2/Modules/zlib/example.o \
./external/python/Python-2.7.2/Modules/zlib/gzio.o \
./external/python/Python-2.7.2/Modules/zlib/infback.o \
./external/python/Python-2.7.2/Modules/zlib/inffast.o \
./external/python/Python-2.7.2/Modules/zlib/inflate.o \
./external/python/Python-2.7.2/Modules/zlib/inftrees.o \
./external/python/Python-2.7.2/Modules/zlib/minigzip.o \
./external/python/Python-2.7.2/Modules/zlib/trees.o \
./external/python/Python-2.7.2/Modules/zlib/uncompr.o \
./external/python/Python-2.7.2/Modules/zlib/zutil.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/zlib/adler32.d \
./external/python/Python-2.7.2/Modules/zlib/compress.d \
./external/python/Python-2.7.2/Modules/zlib/crc32.d \
./external/python/Python-2.7.2/Modules/zlib/deflate.d \
./external/python/Python-2.7.2/Modules/zlib/example.d \
./external/python/Python-2.7.2/Modules/zlib/gzio.d \
./external/python/Python-2.7.2/Modules/zlib/infback.d \
./external/python/Python-2.7.2/Modules/zlib/inffast.d \
./external/python/Python-2.7.2/Modules/zlib/inflate.d \
./external/python/Python-2.7.2/Modules/zlib/inftrees.d \
./external/python/Python-2.7.2/Modules/zlib/minigzip.d \
./external/python/Python-2.7.2/Modules/zlib/trees.d \
./external/python/Python-2.7.2/Modules/zlib/uncompr.d \
./external/python/Python-2.7.2/Modules/zlib/zutil.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/zlib/%.o: ../external/python/Python-2.7.2/Modules/zlib/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


