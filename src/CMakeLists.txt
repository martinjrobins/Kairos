file(GLOB SRC_FILES "*.cpp")
file(GLOB HEADER_FILES "*.h")

if(NOT HAVE_VTK)
	add_definitions(-DHAVE_VTK)
	find_package(VTK REQUIRED)
endif()

FIND_PACKAGE(Boost REQUIRED)

include_directories(${VTK_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++0x")

if(OPTION_STATIC)
	add_library(nsv STATIC ${SRC_FILES} ${HEADER_FILES})
else()
	add_library(nsv SHARED ${SRC_FILES} ${HEADER_FILES})
endif(OPTION_STATIC)

TARGET_LINK_LIBRARIES(nsv ${VTK_LIBRARIES})

