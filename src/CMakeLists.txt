# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.

add_executable (NCsbn NCsbn.cxx)
add_executable (osc_NCsbn osc_NCsbn.cxx)

target_link_libraries (NCsbn SBNfitlib tinyxmllib ${ROOT_LIBRARIES})
target_link_libraries (osc_NCsbn SBNfitlib tinyxmllib ${ROOT_LIBRARIES} ${GSL_LIBRARIES})


