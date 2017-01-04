################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/DrEffect/Animation.cpp \
../src/DrEffect/DrDefinition.cpp \
../src/DrEffect/DrEffect.cpp \
../src/DrEffect/DrFinestra.cpp \
../src/DrEffect/DrOpenGL.cpp \
../src/DrEffect/DrScript.cpp 

OBJS += \
./src/DrEffect/Animation.o \
./src/DrEffect/DrDefinition.o \
./src/DrEffect/DrEffect.o \
./src/DrEffect/DrFinestra.o \
./src/DrEffect/DrOpenGL.o \
./src/DrEffect/DrScript.o 

CPP_DEPS += \
./src/DrEffect/Animation.d \
./src/DrEffect/DrDefinition.d \
./src/DrEffect/DrEffect.d \
./src/DrEffect/DrFinestra.d \
./src/DrEffect/DrOpenGL.d \
./src/DrEffect/DrScript.d 


# Each subdirectory must supply rules for building sources it contributes
src/DrEffect/%.o: ../src/DrEffect/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


