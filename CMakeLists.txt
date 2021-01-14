# -*- cmake -*-
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(HighLANDPD)

# ---- create folders ----------
file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/include")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/dict")

#-------------------------------------------------------
# ----- Options for reading LArSoft files --------------
#-------------------------------------------------------

set (LARSOFT_PRODUCTION PDSPProd2)
set (LARSOFT_DATAMC_DIR mc)

set (LARSOFTREADER_DIR "${CMAKE_SOURCE_DIR}/src/LArSoftReader/src/${LARSOFT_PRODUCTION}")

#-------------------------------------------------------
# ----- Set cmake options ------------------------------
#-------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-g -pipe -Wall -ansi -Wnon-virtual-dtor -fpic -Woverloaded-virtual -std=c++17")

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")

message("CMAKE source dir: ${CMAKE_SOURCE_DIR}")
message("CMAKE binary dir: ${CMAKE_BINARY_DIR}")
message("CMAKE compiler:   ${CMAKE_CXX_COMPILER}")
message("CMAKE_CXX_FLAGS:  ${CMAKE_CXX_FLAGS}")
message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

#-------------------------------------------------------
#------- External packages -----------------------------
#-------------------------------------------------------

set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# Look for the required libraries

#------ ROOT -------------------
find_package(ROOT REQUIRED COMPONENTS Geom EG)
#find_package(ROOT REQUIRED COMPONENTS Geom EG TMVA)
include(${ROOT_USE_FILE})

message("ROOT_USE_FILE:            ${ROOT_USE_FILE}")
message("ROOT include directories: ${ROOT_INCLUDE_DIR}")
message("ROOT libraries:           ${ROOT_LIBRARIES}")

#------ HighLAND -------------------
find_package(HIGHLAND REQUIRED)

message("HighLAND_USE_FILE:            ${HIGHLAND_USE_FILE}")
message("HighLAND include directories: ${HIGHLAND_INCLUDE_DIRS}")
message("HighLAND libraries:           ${HIGHLAND_LIBRARIES}")

#------ TENSORFLOW------------------
find_package(TENSORFLOW REQUIRED)

message("TENSORFLOW include directories: ${TENSORFLOW_INCLUDE_DIRS}")
message("TENSORFLOW libraries:           ${TENSORFLOW_LIBRARIES}")

#------ PROTOBUF------------------
find_package(PROTOBUF REQUIRED)

message("PROTOBUF include directories: ${PROTOBUF_INCLUDE_DIRS}")
message("PROTOBUF libraries:           ${PROTOBUF_LIBRARIES}")

#------ EIGEN------------------
find_package(EIGEN REQUIRED)

message("EIGEN include directories: ${EIGEN_INCLUDE_DIRS}")
message("EIGEN libraries:           ${EIGEN_LIBRARIES}")


#-------------------------------------------------------
#----- include directories -----------------------------
#-------------------------------------------------------

file(GLOB highlandPD_inc ${CMAKE_SOURCE_DIR}/src/*/src
                         ${CMAKE_SOURCE_DIR}/src/*/src/CNN
                         ${CMAKE_SOURCE_DIR}/src/*/src/corrections ${CMAKE_SOURCE_DIR}/src/*/src/systematics
#                         ${LARSOFTREADER_DIR} ${LARSOFTREADER_DIR}/${LARSOFT_DATAMC_DIR}
                         ${CMAKE_BINARY_DIR}/dict)

include_directories(${highlandPD_inc})
message("Include directories: ${INCLUDE_DIRECTORIES}")

#-------------------------------------------------------
#------- ROOT Dictionaries -----------------------------
#-------------------------------------------------------

file(GLOB dict_src ${CMAKE_SOURCE_DIR}/src/*/src/*_LinkDef.h)

foreach(file0 ${dict_src})
  string(FIND ${file0} "/"          pos_ini REVERSE)
  string(SUBSTRING ${file0} 0      ${pos_ini}  dir)
  string(REPLACE "${dir}/"     "" file2 ${file0})
  string(REPLACE _LinkDef.h    "" file3 ${file2})  

  root_generate_dictionary(dict/${file3}_Dict NODEPHEADERS ${file3}.hxx LINKDEF ${file0}  MODULE ${file3} DEPENDENCIES ${ROOT_LIBRARIES})
#  root_generate_dictionary(dict/${file3}_Dict ${dir}/${file3}.hxx LINKDEF ${file0}  MODULE ${file3} DEPENDENCIES ${ROOT_LIBRARIES})

  list(APPEND dict_list ${dir}/${file3}.cxx dict/${file3}_Dict.cxx)
  list(APPEND dict_src  ${CMAKE_BINARY_DIR}/dict/${file3}_Dict.cxx)

  INSTALL(FILES "${CMAKE_BINARY_DIR}/lib${file3}_rdict.pcm" DESTINATION ${CMAKE_SOURCE_DIR}/lib)
  INSTALL(FILES "${CMAKE_BINARY_DIR}/lib${file3}.rootmap"   DESTINATION ${CMAKE_SOURCE_DIR}/lib)
endforeach(file0)

#  root_generate_dictionary(dict/LArSoftReaderProjectDict
##                           NODEPHEADERS LArSoftReaderProjectHeaders.h
#                           ${LARSOFTREADER_DIR}/${LARSOFT_DATAMC_DIR}/LArSoftReaderProjectHeaders.h
#                           LINKDEF ${LARSOFTREADER_DIR}/${LARSOFT_DATAMC_DIR}/LArSoftReaderLinkDef.h
#                           MODULE LArSoftReader
#                           DEPENDENCIES ${ROOT_LIBRARIES})
#  list(APPEND dict_list ${LARSOFTREADER_DIR}/${LARSOFT_DATAMC_DIR}/LArSoftReaderProjectSource.cxx dict/LArSoftReaderProjectDict.cxx)
#  INSTALL(FILES "${CMAKE_BINARY_DIR}/libLArSoftReader_rdict.pcm" DESTINATION ${CMAKE_SOURCE_DIR}/lib)
#  INSTALL(FILES "${CMAKE_BINARY_DIR}/libLArSoftReader.rootmap"   DESTINATION ${CMAKE_SOURCE_DIR}/lib)


#-------------------------------------------------------
#-------- the highlandPD library -------------------------
#-------------------------------------------------------

file(GLOB highlandPD_src  ${CMAKE_SOURCE_DIR}/src/*/src/*.cxx
                          ${CMAKE_SOURCE_DIR}/src/*/src/CNN/*.cxx
                          ${CMAKE_SOURCE_DIR}/src/*/src/corrections/*.cxx ${CMAKE_SOURCE_DIR}/src/*/src/systematics/*.cxx) 
#                          ${LARSOFTREADER_DIR}/LArSoft*.cxx
#                          ${LARSOFTREADER_DIR}/${LARSOFT_DATAMC_DIR}/LArSoftReaderProjectSource.cxx)    

add_library(highlandPD SHARED ${dict_src} ${highlandPD_src} ) 
target_include_directories(highlandPD PRIVATE ${HIGHLAND_INCLUDE_DIRS} ${highlandPD_inc} ${TENSORFLOW_INCLUDE_DIRS} ${PROTOBUF_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS})


target_link_libraries(highlandPD ${ROOT_LIBRARIES} ${HIGHLAND_LIBRARIES})
INSTALL(TARGETS highlandPD DESTINATION ${CMAKE_SOURCE_DIR}/lib)

set(HIGHLAND_PD_LIBRARIES "-L${CMAKE_BINARY_DIR} -lhighlandPD")

#-------------------------------------------------------
#----- highlandPD includes -------------------------------
#-------------------------------------------------------

file(GLOB highlandPD_inc2 ${CMAKE_SOURCE_DIR}/src/*/src/*.hxx
                          ${CMAKE_SOURCE_DIR}/src/*/src/CNN/*.hxx
                          ${CMAKE_SOURCE_DIR}/src/*/src/corrections/*.hxx ${CMAKE_SOURCE_DIR}/src/*/src/systematics/*.hxx) 
#                          ${LARSOFTREADER_DIR}/*.hxx
#                          ${LARSOFTREADER_DIR}/${LARSOFT_DATAMC_DIR}/*.h)

foreach(file0 ${highlandPD_inc2})
   INSTALL(FILES ${file0} DESTINATION ${CMAKE_SOURCE_DIR}/include)
endforeach(file0)

#-------------------------------------------------------
#-------------- executable -----------------------------
#-------------------------------------------------------

#add_executable(CreateMiniTree ${CMAKE_SOURCE_DIR}/src/highland2/highlandIO/app/RunCreateMiniTree.cxx)
#target_include_directories(CreateMiniTree PRIVATE ${highlandPD_inc})
#target_link_libraries(CreateMiniTree ${ROOT_LIBRARIES} "-L${CMAKE_BINARY_DIR} -lhighlandPD")
#INSTALL(TARGETS CreateMiniTree DESTINATION ${CMAKE_SOURCE_DIR}/bin)

add_executable(pion.exe ${CMAKE_SOURCE_DIR}/src/pionAnalysis/app/RunPionAnalysis.cxx)
target_include_directories(pion.exe PRIVATE ${HIGHLAND_INCLUDE_DIRS} ${highlandPD_inc})
target_link_libraries(pion.exe ${ROOT_LIBRARIES} ${HIGHLAND_LIBRARIES} ${HIGHLAND_PD_LIBRARIES})
INSTALL(TARGETS pion.exe DESTINATION ${CMAKE_SOURCE_DIR}/bin)

add_executable(CreateMiniTreePionAna.exe ${CMAKE_SOURCE_DIR}/src/pionAnalysis/app/RunCreateMiniTree.cxx)
target_include_directories(CreateMiniTreePionAna.exe PRIVATE ${HIGHLAND_INCLUDE_DIRS} ${highlandPD_inc})
target_link_libraries(CreateMiniTreePionAna.exe ${ROOT_LIBRARIES} ${HIGHLAND_LIBRARIES} ${HIGHLAND_PD_LIBRARIES})
INSTALL(TARGETS CreateMiniTreePionAna.exe DESTINATION ${CMAKE_SOURCE_DIR}/bin)

add_executable(proton.exe ${CMAKE_SOURCE_DIR}/src/stoppingProtonAnalysis/app/RunStoppingProtonAnalysis.cxx)
target_include_directories(proton.exe PRIVATE ${HIGHLAND_INCLUDE_DIRS} ${highlandPD_inc})
target_link_libraries(proton.exe ${ROOT_LIBRARIES} ${HIGHLAND_LIBRARIES} ${HIGHLAND_PD_LIBRARIES})
INSTALL(TARGETS proton.exe DESTINATION ${CMAKE_SOURCE_DIR}/bin)

add_executable(dedx.exe ${CMAKE_SOURCE_DIR}/src/pdStudies/app/RundEdxStudies.cxx)
target_include_directories(dedx.exe PRIVATE ${HIGHLAND_INCLUDE_DIRS} ${highlandPD_inc})
#target_link_libraries(dedx.exe ${ROOT_LIBRARIES} ${HIGHLAND_LIBRARIES} ${HIGHLAND_PD_LIBRARIES} ${TENSORFLOW_LIBRARIES} ${PROTOBUF_LIBRARIES})
target_link_libraries(dedx.exe ${ROOT_LIBRARIES} ${HIGHLAND_LIBRARIES} ${HIGHLAND_PD_LIBRARIES})
INSTALL(TARGETS dedx.exe DESTINATION ${CMAKE_SOURCE_DIR}/bin)


#-------------------------------------------------------
#------- custom targets --------------------------------
#-------------------------------------------------------

# target to create dictionaries before libraries
add_custom_target(highlandPD_dicts ALL SOURCES ${dict_list})

# target to change the library name in the rootmap files
add_custom_target(ch_rootmap COMMAND ${CMAKE_SOURCE_DIR}/scripts/chrootmap ${CMAKE_BINARY_DIR})

# target to install header files in the include folder
add_custom_target(include ALL SOURCES ${highlandPD_inc2} ${CMAKE_SOURCE_DIR})

#-------------------------------------------------------
#------- dependencies ----------------------------------
#-------------------------------------------------------

add_dependencies(dedx.exe CreateMiniTreePionAna.exe )
add_dependencies(CreateMiniTreePionAna.exe proton.exe)
add_dependencies(proton.exe                pion.exe)
add_dependencies(pion.exe                  highlandPD)
add_dependencies(highlandPD                ch_rootmap)
add_dependencies(ch_rootmap                highlandPD_dicts)
add_dependencies(highlandPD_dicts          include)

#-------------------------------------------------------
#------- other -----------------------------------------
#-------------------------------------------------------

#---Configure and install various files neded later and for clients -----------------------------
include(HighlandPDConfiguration)
