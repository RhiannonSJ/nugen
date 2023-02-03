#art_make(MODULE_LIBRARIES nusimdata_SimulationBase
#                          larcoreobj_SummaryData
#                          art_Framework_Services_Optional_RandomNumberGenerator_service
#                          art_Framework_Services_Registry
#                          ${MF_MESSAGELOGGER}
#                          ${MF_UTILITIES}
#                          ${FHICLCPP}
#                          ${CETLIB}
#                          ${LOG4CPP}
#                          ${XML2}
#                          ${PYTHIA6}
#                          ${LHAPDF}
#                          ${GENIE_LIB_LIST}
#                          ${CLHEP}
#                          ${DK2NU_TREE}
#                          ${ROOT_BASIC_LIB_LIST}
#                          ${ROOT_EG}
#                          ${ROOT_TREEPLAYER}
#                          ${ROOT_FFTW}
#                          ${ROOT_REFLEX}
#                          ${ROOTSYS}/lib/libEGPythia6.so
#                          ${ROOT_GUI}
#                          ${ROOT_EGPythia6}
#                          nugen_EventGeneratorBase_LArModules
#                        )
cet_build_plugin(AddGenieEventsToLAr art::module
                 LIBRARIES REG nusimdata::SimulationBase
                        nugen::EventGeneratorBase_GENIE
                        larcoreobj_SummaryData
                        art::Framework_Services_Optional_RandomNumberGenerator_service
                        art::Framework_Services_Registry
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        log4cpp::log4cpp
                        LibXml2::LibXml2
                        Pythia6::Pythia6
                        LHAPDF::LHAPDF
                        ${GENIE_LIB_LIST}
                        GSL::gsl  # FWNUM relies on GSL
                        BLAS::BLAS
                        ROOT::EGPythia6
                        ${ROOT_EVE_LIB_LIST}
                        CLHEP::CLHEP
                        dk2nu::Tree
                        ROOT::Core
                        ROOT::MathMore
                        ROOT::EG
                        ROOT::TreePlayer
                        ROOT::FFTW
                        ROOT::EGPythia6
                        ROOT::Gui
                        )

install_headers()
install_fhicl()
install_source( EXTRAS README )
