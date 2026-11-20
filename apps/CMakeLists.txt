
add_executable(kff-tools
	kfftools.cpp
)

target_link_libraries(kff-tools
	kfftl
)

# Version of kff-tools, given by --version
target_compile_definitions(kff-tools PRIVATE KFF_TOOLS_VERSION="${kff-tools_VERSION}")
