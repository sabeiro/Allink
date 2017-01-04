################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/VarData/Cubo.o \
../src/VarData/UsaMatematica.o \
../src/VarData/VarDatFile.o \
../src/VarData/VarData.o \
../src/VarData/VarDataBackFold.o \
../src/VarData/VarDataCGAL.o \
../src/VarData/VarDataComm.o \
../src/VarData/VarDataCreate.o \
../src/VarData/VarDataEl.o \
../src/VarData/VarDataExp.o \
../src/VarData/VarDataInterp.o \
../src/VarData/VarDataMarchCubes.o \
../src/VarData/VarDataProfile.o \
../src/VarData/VarDataRead.o \
../src/VarData/VarDataString.o \
../src/VarData/VarDataWrite.o 

CPP_SRCS += \
../src/VarData/Cubo.cpp \
../src/VarData/Polymers.cpp \
../src/VarData/UsaMatematica.cpp \
../src/VarData/VarDatFile.cpp \
../src/VarData/VarData.cpp \
../src/VarData/VarDataBackFold.cpp \
../src/VarData/VarDataCGAL.cpp \
../src/VarData/VarDataComm.cpp \
../src/VarData/VarDataContour.cpp \
../src/VarData/VarDataCreate.cpp \
../src/VarData/VarDataEl.cpp \
../src/VarData/VarDataExp.cpp \
../src/VarData/VarDataInterp.cpp \
../src/VarData/VarDataMarchCubes.cpp \
../src/VarData/VarDataPos.cpp \
../src/VarData/VarDataProfile.cpp \
../src/VarData/VarDataRead.cpp \
../src/VarData/VarDataString.cpp \
../src/VarData/VarDataWrite.cpp 

OBJS += \
./src/VarData/Cubo.o \
./src/VarData/Polymers.o \
./src/VarData/UsaMatematica.o \
./src/VarData/VarDatFile.o \
./src/VarData/VarData.o \
./src/VarData/VarDataBackFold.o \
./src/VarData/VarDataCGAL.o \
./src/VarData/VarDataComm.o \
./src/VarData/VarDataContour.o \
./src/VarData/VarDataCreate.o \
./src/VarData/VarDataEl.o \
./src/VarData/VarDataExp.o \
./src/VarData/VarDataInterp.o \
./src/VarData/VarDataMarchCubes.o \
./src/VarData/VarDataPos.o \
./src/VarData/VarDataProfile.o \
./src/VarData/VarDataRead.o \
./src/VarData/VarDataString.o \
./src/VarData/VarDataWrite.o 

CPP_DEPS += \
./src/VarData/Cubo.d \
./src/VarData/Polymers.d \
./src/VarData/UsaMatematica.d \
./src/VarData/VarDatFile.d \
./src/VarData/VarData.d \
./src/VarData/VarDataBackFold.d \
./src/VarData/VarDataCGAL.d \
./src/VarData/VarDataComm.d \
./src/VarData/VarDataContour.d \
./src/VarData/VarDataCreate.d \
./src/VarData/VarDataEl.d \
./src/VarData/VarDataExp.d \
./src/VarData/VarDataInterp.d \
./src/VarData/VarDataMarchCubes.d \
./src/VarData/VarDataPos.d \
./src/VarData/VarDataProfile.d \
./src/VarData/VarDataRead.d \
./src/VarData/VarDataString.d \
./src/VarData/VarDataWrite.d 


# Each subdirectory must supply rules for building sources it contributes
src/VarData/%.o: ../src/VarData/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


