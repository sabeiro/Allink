################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/Visualizza/ElementiGrafici.o \
../src/Visualizza/ElementiGraficiComm.o \
../src/Visualizza/ElementiGraficiDis.o \
../src/Visualizza/ElementiGraficiEl.o \
../src/Visualizza/Finestra.o \
../src/Visualizza/VarDatFile.o \
../src/Visualizza/Visualizza.o \
../src/Visualizza/moc_ElementiGrafici.o 

CPP_SRCS += \
../src/Visualizza/ElementiGrafici.cpp \
../src/Visualizza/ElementiGraficiComm.cpp \
../src/Visualizza/ElementiGraficiDis.cpp \
../src/Visualizza/ElementiGraficiEl.cpp \
../src/Visualizza/Finestra.cpp \
../src/Visualizza/Visualizza.cpp \
../src/Visualizza/main.cpp \
../src/Visualizza/moc_ElementiGrafici.cpp 

OBJS += \
./src/Visualizza/ElementiGrafici.o \
./src/Visualizza/ElementiGraficiComm.o \
./src/Visualizza/ElementiGraficiDis.o \
./src/Visualizza/ElementiGraficiEl.o \
./src/Visualizza/Finestra.o \
./src/Visualizza/Visualizza.o \
./src/Visualizza/main.o \
./src/Visualizza/moc_ElementiGrafici.o 

CPP_DEPS += \
./src/Visualizza/ElementiGrafici.d \
./src/Visualizza/ElementiGraficiComm.d \
./src/Visualizza/ElementiGraficiDis.d \
./src/Visualizza/ElementiGraficiEl.d \
./src/Visualizza/Finestra.d \
./src/Visualizza/Visualizza.d \
./src/Visualizza/main.d \
./src/Visualizza/moc_ElementiGrafici.d 


# Each subdirectory must supply rules for building sources it contributes
src/Visualizza/%.o: ../src/Visualizza/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


