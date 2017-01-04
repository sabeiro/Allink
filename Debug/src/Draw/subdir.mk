################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/Draw/Draw.o \
../src/Draw/DrawControl.o \
../src/Draw/DrawDefinition.o \
../src/Draw/DrawFile.o \
../src/Draw/DrawScene.o \
../src/Draw/ProvaDraw.o 

CPP_SRCS += \
../src/Draw/DrImage.cpp \
../src/Draw/Draw.cpp \
../src/Draw/DrawControl.cpp \
../src/Draw/DrawDefinition.cpp \
../src/Draw/DrawFile.cpp \
../src/Draw/DrawScene.cpp \
../src/Draw/ProvaDraw.cpp \
../src/Draw/ReadWritePng.cpp 

OBJS += \
./src/Draw/DrImage.o \
./src/Draw/Draw.o \
./src/Draw/DrawControl.o \
./src/Draw/DrawDefinition.o \
./src/Draw/DrawFile.o \
./src/Draw/DrawScene.o \
./src/Draw/ProvaDraw.o \
./src/Draw/ReadWritePng.o 

CPP_DEPS += \
./src/Draw/DrImage.d \
./src/Draw/Draw.d \
./src/Draw/DrawControl.d \
./src/Draw/DrawDefinition.d \
./src/Draw/DrawFile.d \
./src/Draw/DrawScene.d \
./src/Draw/ProvaDraw.d \
./src/Draw/ReadWritePng.d 


# Each subdirectory must supply rules for building sources it contributes
src/Draw/%.o: ../src/Draw/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


