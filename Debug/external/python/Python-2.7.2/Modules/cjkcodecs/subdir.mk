################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_cn.c \
../external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_hk.c \
../external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_iso2022.c \
../external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_jp.c \
../external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_kr.c \
../external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_tw.c \
../external/python/Python-2.7.2/Modules/cjkcodecs/multibytecodec.c 

OBJS += \
./external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_cn.o \
./external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_hk.o \
./external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_iso2022.o \
./external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_jp.o \
./external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_kr.o \
./external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_tw.o \
./external/python/Python-2.7.2/Modules/cjkcodecs/multibytecodec.o 

C_DEPS += \
./external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_cn.d \
./external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_hk.d \
./external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_iso2022.d \
./external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_jp.d \
./external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_kr.d \
./external/python/Python-2.7.2/Modules/cjkcodecs/_codecs_tw.d \
./external/python/Python-2.7.2/Modules/cjkcodecs/multibytecodec.d 


# Each subdirectory must supply rules for building sources it contributes
external/python/Python-2.7.2/Modules/cjkcodecs/%.o: ../external/python/Python-2.7.2/Modules/cjkcodecs/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


