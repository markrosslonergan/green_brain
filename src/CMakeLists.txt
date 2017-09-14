# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.

add_executable (NCsbn NCsbn.cxx)

target_link_libraries (NCsbn SBNfitlib tinyxmllib ${ROOT_LIBRARIES})

