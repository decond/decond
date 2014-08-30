cmake_minimum_required (VERSION 2.6)
project (econd)
enable_language (Fortran)

set (VERSION_MAJOR "1")
set (VERSION_MINOR "0")
set (VERSION_PATCH "0")
set (VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

add_subdirectory (src)

# install executables and scripts
#set (EXECUTABLES "NORMA.exe" ${NMPROGRAMS})
#set (SCRIPTS "gen_pert.sh" "pert_multi_mode.sh")
#install (TARGETS ${EXECUTABLES} RUNTIME DESTINATION "bin")
#install (PROGRAMS ${SCRIPTS} DESTINATION "bin") 

