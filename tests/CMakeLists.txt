add_executable(FinalProjectTests
        "your_tests.cpp"
        "bvh_tests.cpp"
        "interpolation_tests.cpp"
        "texture_tests.cpp"
        "light_tests.cpp"
        "depth_of_field_tests.cpp"
)

target_compile_features(FinalProjectTests PRIVATE cxx_std_20)
set_project_warnings(FinalProjectTests)
target_link_libraries(FinalProjectTests PRIVATE CGFramework FinalProjectLib Catch2WithMain)
