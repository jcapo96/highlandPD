# -*- cmake -*-
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(pionAnalysis)

# ---- create folders ----------
file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/include")

#-------------------------------------------------------
# ----- Set cmake options ------------------------------
#-------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-g -pipe -Wall -ansi -Wnon-virtual-dtor -fpic -Woverloaded-virtual -std=c++17")

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")

#-------------------------------------------------------
#------- External packages -----------------------------
#-------------------------------------------------------

set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# Look for the required libraries

find_package(ROOT          REQUIRED COMPONENTS Geom EG)
find_package(HIGHLAND      REQUIRED)
find_package(HIGHLAND_PDSP REQUIRED)

include(${ROOT_USE_FILE})

include(cmake/scripts/CMakeLists.txt)

#-------------------------------------------------------
#------- ROOT Dictionaries -----------------------------
#-------------------------------------------------------

#root_generate_dictionary(PionAnaDataClasses_Dict
#                        NODEPHEADERS PionAnaDataClasses.hxx
#                         ${CMAKE_SOURCE_DIR}/src/PionAnaDataClasses.hxx
#                         LINKDEF ${CMAKE_SOURCE_DIR}/src/PionAnaDataClasses_LinkDef.h MODULE PionAnaDataClasses DEPENDENCIES ${ROOT_LIBRARIES})

# target to create dictionaries before libraries
#add_custom_target(dicts ALL SOURCES PionAnaDataClasses_Dict.cxx)

#-------------------------------------------------------
#-------- the pionAnalysis library -------------------------
#-------------------------------------------------------

file(GLOB pionana_src  ${CMAKE_SOURCE_DIR}/src/*.cxx
                       ${CMAKE_SOURCE_DIR}/src/corrections/*.cxx 
                       ${CMAKE_SOURCE_DIR}/src/systematics/*.cxx)
#                       ${CMAKE_BINARY_DIR}/PionAnaDataClasses_Dict.cxx)

add_library(pionana SHARED ${pionana_src}) 
target_include_directories(pionana PRIVATE ${HIGHLAND_INCLUDE_DIRS} ${HIGHLAND_PDSP_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/corrections ${CMAKE_SOURCE_DIR}/src/systematics)
target_link_libraries(pionana ${ROOT_LIBRARIES} ${HIGHLAND_LIBRARIES} ${HIGHLAND_PDSP_LIBRARIES})
INSTALL(TARGETS pionana DESTINATION ${CMAKE_SOURCE_DIR}/lib)

set(PIONANA_LIBRARIES "-L${CMAKE_BINARY_DIR} -lpionana")

#-------------------------------------------------------
#-------------- executables ----------------------------
#-------------------------------------------------------

add_executable(pion.exe ${CMAKE_SOURCE_DIR}/app/RunPionAnalysis.cxx)
target_include_directories(pion.exe PRIVATE ${HIGHLAND_INCLUDE_DIRS} ${HIGHLAND_PDSP_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(pion.exe ${ROOT_LIBRARIES} ${HIGHLAND_LIBRARIES} ${HIGHLAND_PDSP_LIBRARIES} ${PIONANA_LIBRARIES})
INSTALL(TARGETS pion.exe DESTINATION ${CMAKE_SOURCE_DIR}/bin)

add_executable(CreateMiniTreePionAna.exe ${CMAKE_SOURCE_DIR}/app/RunCreateMiniTree.cxx)
target_include_directories(CreateMiniTreePionAna.exe PRIVATE ${HIGHLAND_INCLUDE_DIRS} ${HIGHLAND_PDSP_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(CreateMiniTreePionAna.exe ${ROOT_LIBRARIES} ${HIGHLAND_LIBRARIES} ${HIGHLAND_PDSP_LIBRARIES} ${PIONANA_LIBRARIES})
INSTALL(TARGETS CreateMiniTreePionAna.exe DESTINATION ${CMAKE_SOURCE_DIR}/bin)
 
#-------------------------------------------------------
#------- dependencies ----------------------------------
#-------------------------------------------------------

add_dependencies(CreateMiniTreePionAna.exe pion.exe)
add_dependencies(pion.exe                  pionana)
#add_dependencies(pionana                   dicts)
