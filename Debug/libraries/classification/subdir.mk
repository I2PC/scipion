################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../libraries/classification/analyze_cluster.cpp \
../libraries/classification/batch_som.cpp \
../libraries/classification/code_book.cpp \
../libraries/classification/fcmeans.cpp \
../libraries/classification/fkcn.cpp \
../libraries/classification/fuzzy_code_book.cpp \
../libraries/classification/fuzzy_som.cpp \
../libraries/classification/gaussian_kerdensom.cpp \
../libraries/classification/kSVD.cpp \
../libraries/classification/kerdensom.cpp \
../libraries/classification/knn_classifier.cpp \
../libraries/classification/map.cpp \
../libraries/classification/naive_bayes.cpp \
../libraries/classification/pca.cpp \
../libraries/classification/sammon.cpp \
../libraries/classification/som.cpp \
../libraries/classification/svm.cpp \
../libraries/classification/svm_classifier.cpp \
../libraries/classification/training_vector.cpp 

OBJS += \
./libraries/classification/analyze_cluster.o \
./libraries/classification/batch_som.o \
./libraries/classification/code_book.o \
./libraries/classification/fcmeans.o \
./libraries/classification/fkcn.o \
./libraries/classification/fuzzy_code_book.o \
./libraries/classification/fuzzy_som.o \
./libraries/classification/gaussian_kerdensom.o \
./libraries/classification/kSVD.o \
./libraries/classification/kerdensom.o \
./libraries/classification/knn_classifier.o \
./libraries/classification/map.o \
./libraries/classification/naive_bayes.o \
./libraries/classification/pca.o \
./libraries/classification/sammon.o \
./libraries/classification/som.o \
./libraries/classification/svm.o \
./libraries/classification/svm_classifier.o \
./libraries/classification/training_vector.o 

CPP_DEPS += \
./libraries/classification/analyze_cluster.d \
./libraries/classification/batch_som.d \
./libraries/classification/code_book.d \
./libraries/classification/fcmeans.d \
./libraries/classification/fkcn.d \
./libraries/classification/fuzzy_code_book.d \
./libraries/classification/fuzzy_som.d \
./libraries/classification/gaussian_kerdensom.d \
./libraries/classification/kSVD.d \
./libraries/classification/kerdensom.d \
./libraries/classification/knn_classifier.d \
./libraries/classification/map.d \
./libraries/classification/naive_bayes.d \
./libraries/classification/pca.d \
./libraries/classification/sammon.d \
./libraries/classification/som.d \
./libraries/classification/svm.d \
./libraries/classification/svm_classifier.d \
./libraries/classification/training_vector.d 


# Each subdirectory must supply rules for building sources it contributes
libraries/classification/%.o: ../libraries/classification/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


