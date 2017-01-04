################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/Matematica/Matematica.o \
../src/Matematica/MatematicaAlgebra.o \
../src/Matematica/MatematicaChar.o \
../src/Matematica/MatematicaFilter.o \
../src/Matematica/MatematicaFunc.o \
../src/Matematica/MatematicaInterp.o \
../src/Matematica/MatematicaQuaternion.o \
../src/Matematica/MatematicaSign.o \
../src/Matematica/MatematicaVect.o \
../src/Matematica/mt19937ar.o \
../src/Matematica/randomlib.o \
../src/Matematica/rnorrexp.o 

CPP_SRCS += \
../src/Matematica/Matematica.cpp \
../src/Matematica/MatematicaAlgebra.cpp \
../src/Matematica/MatematicaChar.cpp \
../src/Matematica/MatematicaFilter.cpp \
../src/Matematica/MatematicaFunc.cpp \
../src/Matematica/MatematicaInterp.cpp \
../src/Matematica/MatematicaQuaternion.cpp \
../src/Matematica/MatematicaSign.cpp \
../src/Matematica/MatematicaVect.cpp \
../src/Matematica/mt19937ar.cpp 

C_SRCS += \
../src/Matematica/Marsenne.c \
../src/Matematica/randomlib.c \
../src/Matematica/rnorrexp.c 

OBJS += \
./src/Matematica/Marsenne.o \
./src/Matematica/Matematica.o \
./src/Matematica/MatematicaAlgebra.o \
./src/Matematica/MatematicaChar.o \
./src/Matematica/MatematicaFilter.o \
./src/Matematica/MatematicaFunc.o \
./src/Matematica/MatematicaInterp.o \
./src/Matematica/MatematicaQuaternion.o \
./src/Matematica/MatematicaSign.o \
./src/Matematica/MatematicaVect.o \
./src/Matematica/mt19937ar.o \
./src/Matematica/randomlib.o \
./src/Matematica/rnorrexp.o 

C_DEPS += \
./src/Matematica/Marsenne.d \
./src/Matematica/randomlib.d \
./src/Matematica/rnorrexp.d 

CPP_DEPS += \
./src/Matematica/Matematica.d \
./src/Matematica/MatematicaAlgebra.d \
./src/Matematica/MatematicaChar.d \
./src/Matematica/MatematicaFilter.d \
./src/Matematica/MatematicaFunc.d \
./src/Matematica/MatematicaInterp.d \
./src/Matematica/MatematicaQuaternion.d \
./src/Matematica/MatematicaSign.d \
./src/Matematica/MatematicaVect.d \
./src/Matematica/mt19937ar.d 


# Each subdirectory must supply rules for building sources it contributes
src/Matematica/%.o: ../src/Matematica/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/Matematica/%.o: ../src/Matematica/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


