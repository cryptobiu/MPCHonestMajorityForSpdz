cmake_minimum_required(VERSION 2.8)
project(MPCHonestMajorityForSpdz)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g  -std=c++11 -maes -mpclmul -mbmi2 -fPIC" )



INCLUDE_DIRECTORIES($ENV{HOME}/ntl-11.1.0/include $ENV{HOME} 
        $ENV{HOME}/libscapi/lib/OTExtensionBristol $ENV{HOME}/libscapi/install/include/)

set(SOURCE_FILES HIM.cpp VDM.cpp MPCCommunication.cpp MPCCommunication.h
        Protocol.cpp Protocol.h TemplateField.h TemplateField.cpp ProtocolTimer.cpp ProtocolTimer.h
        HonestMultAbstract.h GRRHonestMult.h DNHonestMult.h
        ZpMersenneLongElement.cpp ZpMersenneLongElement.h Mersenne127.cpp Mersenne127.h)


link_directories(/usr/ssl/lib/ $ENV{HOME}/libscapi/install/lib )

#add_executable(MPCHonestMajorityForSpdz ${SOURCE_FILES})
add_library(MPCHonestMajorityForSpdz ${SOURCE_FILES})


TARGET_LINK_LIBRARIES(MPCHonestMajorityForSpdz $ENV{HOME}/ntl-11.1.0/src/ntl.a
        $ENV{HOME}/libscapi/libscapi.a -lgmp gmpxx
        OTExtensionBristol
        $ENV{HOME}/libscapi/install/lib/libsimpleot.a
        libboost_system.a libboost_thread.a pthread crypto dl ssl z )


