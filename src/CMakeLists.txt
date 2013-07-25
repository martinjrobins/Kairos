file(GLOB Kairos_SOURCES "*.cpp")
file(GLOB Kairos_HEADERS "*.h")
include_directories(${Kairos_INCLUDE_DIRECTORIES})
message(STATUS ${Kairos_INCLUDE_DIRECTORIES})

add_library(Kairos SHARED ${Kairos_SOURCES} "../smoldyn/rxnparam.c")
TARGET_LINK_LIBRARIES(Kairos ${VTK_LIBRARIES})

install (TARGETS Kairos DESTINATION lib)
install (FILES ${Kairos_HEADERS} DESTINATION include)
install (FILES ../smoldyn/rxnparam.h DESTINATION include)
install (DIRECTORY ${Kairos_SOURCE_DIR}/eigen3/ DESTINATION include)

