cet_make_library(LIBRARY_NAME EventGeneratorBase INTERFACE SOURCE evgenbase.h)

install_headers()
install_fhicl()
install_source()

add_subdirectory(GENIE)
add_subdirectory(GiBUU)
add_subdirectory(Modules)
#add_subdirectory(test)
